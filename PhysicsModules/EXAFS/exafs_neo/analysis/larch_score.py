import os
import csv
import numpy as np
from typing import List, Dict, Tuple
import matplotlib.pyplot as plt
# Larch imports
from larch.xafs import autobk, feffdat, xftf
from larch.io import read_ascii
from larch import Interpreter

os.environ["NUMEXPR_MAX_THREADS"] = "16"


def _apply_xftf(group, kmin, kmax, kweight, larch_env):
    """Helper function to apply xftf consistently."""
    xftf(
        group.k,
        group.chi,
        kmin=kmin,
        kmax=kmax,
        dk=4,
        window="hanning",
        kweight=kweight,
        group=group,
        _larch=larch_env,
    )


def larch_init(csv_sub: str, params: Dict) -> Tuple[np.ndarray, object, Dict, Interpreter]:
    """Initialize larch and process base experimental data."""
    mylarch = Interpreter()

    Kmin = params["Kmin"]
    Kmax = params["Kmax"]
    deltak = params["deltak"]

    SMALL = int(Kmin / deltak)
    BIG = int(Kmax / deltak)
    MID = int(BIG - SMALL + 1)
    RBKG = params["rbkg"]
    KWEIGHT = params["kweight"]
    BKGKW = params["bkgkw"]
    BKGKMAX = params["bkgkmax"]

    base = params["base"]
    csv_path = os.path.join(base, csv_sub)

    g = read_ascii(csv_path)
    best = read_ascii(csv_path)
    sumgroup = read_ascii(csv_path)

    if not hasattr(g, 'chi'):
        autobk(g, rbkg=RBKG, kweight=BKGKW, kmax=BKGKMAX, _larch=mylarch)
        autobk(best, rbkg=RBKG, _larch=mylarch)
        autobk(sumgroup, rbkg=RBKG, _larch=mylarch)

    # Convert to integer array for fast indexing later
    intervalK = np.linspace(SMALL, BIG, MID).astype(int)

    # Apply xftf to all groups
    for group in (g, best, sumgroup):
        _apply_xftf(group, Kmin, Kmax, KWEIGHT, mylarch)

    # Store state cleanly in params instead of globals
    params.update({
        "SMALL": SMALL,
        "BIG": BIG,
        "intervalK": intervalK,
        "best_group": best
    })

    return g.chi, g, params, mylarch


def flatten_2d_list(array: List) -> List:
    """Flatten a list that may contain mixed items and sublists."""
    flat_list = []
    for item in array:
        if isinstance(item, list):
            flat_list.extend(item)
        else:
            flat_list.append(item)
    return flat_list


def fitness_individual(exp:np.ndarray, arr, full_paths, params, mylarch, g, plot=False, export=False, fig_gui=None):
    """Fitness for individual scattering paths."""
    intervalK = params["intervalK"]
    best = params["best_group"]
    Kmax = params["Kmax"]
    Kweight = params["kweight"]
    SMALL = params["SMALL"]
    BIG = params["BIG"]
    base = params["base"]

    num_comp = len(params["front"])
    compounds_list = params["front"]

    fronts = [os.path.join(base, i) for i in compounds_list] if num_comp > 1 else [
        os.path.join(base, params["front"][0])]

    yTotal = np.zeros(401)
    flat_paths = flatten_2d_list(full_paths)
    export_paths = np.zeros((2 * len(flat_paths), 401))

    ax = None
    if plot:
        fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(7, 6))
    if fig_gui is not None:
        ax = fig_gui.add_subplot(111)

    iterator = 0
    path = None  # To hold the last path object

    for i in range(num_comp):
        paths = full_paths[i] if num_comp > 1 else full_paths
        front = fronts[i] if num_comp > 1 else fronts[0]

        for j, p in enumerate(paths):
            filename = f"{front}{p:04d}.dat"

            path = feffdat.feffpath(
                filename,
                s02=str(arr[j, 0]),
                e0=str(arr[j, 1]),
                sigma2=str(arr[j, 2]),
                deltar=str(arr[j, 3]),
                _larch=mylarch,
            )
            feffdat.path2chi(path, _larch=mylarch)

            # Vectorized addition
            yTotal[intervalK] += path.chi[intervalK]

            if ax is not None:
                offset = 5
                ax.plot(path.k, path.chi * path.k ** 2.0 + offset * (iterator + 1), label="Path")

            if export:
                export_paths[2 * iterator, :] = path.k
                export_paths[2 * iterator + 1, :] = path.chi * path.k ** 2.0

            iterator += 1

    if ax is not None:
        ax.set_xlabel("k ($\AA^{-1}$)")
        ax.set_ylabel("k$^{2}$ ($\chi(k)\AA^{-1}$)")
        ax.set_ylim(-10, len(flat_paths) * 5 + 5)
        ax.set_xlim(0, Kmax + 1)
        ax.plot(g.k, g.chi * g.k ** 2, "r--", label="Data")
        ax.plot(path.k[SMALL:BIG], yTotal[SMALL:BIG] * path.k[SMALL:BIG] ** 2, "b--", label="GA")
        ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")

    best.chi = yTotal
    best.k = path.k
    _apply_xftf(best, params["Kmin"], Kmax, Kweight, mylarch)

    # Vectorized loss calculation (Square Loss)
    loss = np.sum(((yTotal[intervalK] * g.k[intervalK] ** Kweight) - (exp[intervalK] * g.k[intervalK] ** Kweight)) ** 2)

    return path, yTotal.tolist(), best, loss, export_paths


def fitness(exp:np.ndarray, arr, full_paths, params, mylarch, g, verbose=False) -> tuple:
    """Fitness of full scattering paths."""
    base = params["base"]
    intervalK = params["intervalK"]
    best = params["best_group"]
    Kmin = params["Kmin"]
    Kmax = params["Kmax"]
    Kweight = params["kweight"]

    compounds_list = params["front"]
    num_comp = len(compounds_list)
    fronts = [os.path.join(base, i) for i in compounds_list] if num_comp > 1 else [
        os.path.join(base, compounds_list[0])]

    yTotal = np.zeros(401)
    arr_r, og_arr_r = [], []
    array_str = "---------------------\n"

    if verbose:
        print(f"Number of Components: {num_comp}")

    sum_list = [0]
    if num_comp > 1:
        for i in range(num_comp):
            sum_list.append(sum_list[-1] + len(full_paths[i]))

    path = None
    for i in range(num_comp):
        paths = full_paths[i] if num_comp > 1 else full_paths
        front = fronts[i] if num_comp > 1 else fronts[0]

        for j, p in enumerate(paths):
            filename = f"{front}{p:04d}.dat"
            k = sum_list[i] + j

            path = feffdat.feffpath(
                filename,
                s02=str(arr[k][0]),
                e0=str(arr[k][1]),
                sigma2=str(arr[k][2]),
                deltar=str(arr[k][3]),
                _larch=mylarch,
            )
            feffdat.path2chi(path, larch=mylarch)

            if verbose:
                print(
                    f"Path {p} {float(path.s02):.3f} {float(path.e0):.3f} {float(path.sigma2):.3f} {float(path.reff + arr[j, 3]):.3f}")

            arr_r.append([
                float(path.s02), float(path.e0), float(path.sigma2),
                float(path.reff + arr[j, 3]), float(path.degen), float(path.nleg), path.geom
            ])
            og_arr_r.append(path.reff)

            # Vectorized addition
            yTotal[intervalK] += path.chi[intervalK]

    best.chi = yTotal
    best.k = path.k
    _apply_xftf(best, Kmin, Kmax, Kweight, mylarch)

    # Vectorized loss calculation
    loss = np.sum(((yTotal[intervalK] * g.k[intervalK] ** Kweight) - (exp[intervalK] * g.k[intervalK] ** Kweight)) ** 2)

    return path, yTotal.tolist(), best, loss, arr_r, og_arr_r, array_str


def construct_bestfit_mat(gk, gchi, pathk, yTotal, small, big):
    """Construct best fit matrix"""
    exp_data = np.column_stack((gk, gchi * gk ** 2))
    fit_data = np.column_stack((pathk[small:big], np.array(yTotal)[small:big] * pathk[small:big] ** 2))
    return exp_data, fit_data


def write_bestFit_csv(gk, gchi, pathk, yTotal, small, big, name="bestFit.csv", header_base="Sample"):
    """Write out bestFit csv for plotting in Igor Format"""
    exp_data, fit_data = construct_bestfit_mat(gk, gchi, pathk, yTotal, small, big)

    with open(name, mode="w", newline="", encoding="utf-8") as write_file:
        writer = csv.writer(write_file)
        writer.writerow([f"data_{header_base}.k", f"data_{header_base}.chi2"])
        writer.writerows(exp_data)
        writer.writerow([])
        writer.writerow([f"fit_{header_base}.k", f"fit_{header_base}.chi2"])
        writer.writerows(fit_data)


def write_individual_csv(gk, gchi, pathk, yTotal, small, big, export_path, paths, name="Individual.csv",
                         header_base="Sample"):
    """Write individual csv"""
    exp_data, fit_data = construct_bestfit_mat(gk, gchi, pathk, yTotal, small, big)

    with open(name, mode="w", newline="", encoding="utf-8") as write_file:
        writer = csv.writer(write_file)
        writer.writerow([f"data_{header_base}.k", f"data_{header_base}.chi2"])
        writer.writerows(exp_data)
        writer.writerow([])
        writer.writerow([f"fit_{header_base}.k", f"fit_{header_base}.chi2"])
        writer.writerows(fit_data)
        writer.writerow([])

        full_paths = flatten_2d_list(paths)
        for i, p in enumerate(full_paths):
            writer.writerow([f"path_{p}_{header_base}.k", f"path_{p}_{header_base}.chi2"])
            writer.writerows(export_path[(2 * i, 2 * i + 1), :].T)
            writer.writerow([])


def write_result_n_err(full_mat_var, err, name="bestfit_err.csv", header_base="Sample"):
    """Write Result and Error"""
    assert len(full_mat_var) == len(err)
    out_mat = np.column_stack((full_mat_var, err))

    with open(name, mode="w", newline="", encoding="utf-8") as write_file:
        writer = csv.writer(write_file)
        writer.writerow([f"BestFit_{header_base}", f"Err_{header_base}"])
        writer.writerows(out_mat)


def generate_labels(path_list):
    """Generate labels"""
    label, s02_label, sigma2_label, deltaR_label = [], [], [], []

    for i, p in enumerate(path_list):
        s02 = f"s02_{p}"
        sig = f"sigma_{p}"
        dr = f"deltaR_{p}"

        label.append(s02)
        if i == 0: label.append("e0")
        label.extend([sig, dr])

        s02_label.append(s02)
        sigma2_label.append(sig)
        deltaR_label.append(dr)

    return label, s02_label, sigma2_label, deltaR_label


def construct_full_err(err):
    """Construct full err matrix"""
    e0 = err[1]
    err_temp = np.delete(err, 1).reshape((-1, 3))
    e0_arr = np.full((err_temp.shape[0], 1), e0)
    return np.hstack((err_temp[:, :1], e0_arr, err_temp[:, 1:]))


def convert_to_str(val, prec) -> str:
    """convert to string of certain precision"""
    return f"{round(val, prec):.{prec}f}"


def convert_label(select_bestfit_r):
    """Convert labels"""
    return "-".join(item[0] for item in select_bestfit_r)


def cal_err_prec(val_arr):
    """Calculate Error Precision"""
    return [abs(int(np.log10(abs(v))) - 1) if v != 0 else 0 for v in val_arr]


def temp_round_deltaR(val, tol=1e-6) -> float:
    """Round value based on preset tolerance"""
    return 0.001 if val < tol else np.round(val, 3)


def latex_table(paths, best_Fit_r, err_full, chi2r, og_arr_r):
    """convert to latex table output"""
    full_paths = flatten_2d_list(paths)
    label_arr = [convert_label(r[6]) for r in best_Fit_r]
    nleg_arr = [str(int(r[5])) for r in best_Fit_r]

    latex_table_str = f"""
\\begin{{table}}[]
    \\centering
        \\footnotesize
            \\caption{{$\chi_r$ = {chi2r}}}                
            \\begin{{tabular}}{{ccccccccc}}
                \\hline
                \\vspace{{0.05 in}}
                    Path \\# & $N$ & $S_0^2$ & $\\Delta$E$_0$ (eV) & $\\sigma^{{2}}$ (\\AA$^2$) &  R (\\AA) & Fitted R (\\AA)& Legs & Labels\\\\
                    \\hline
"""

    for i, p in enumerate(full_paths):
        prec = cal_err_prec(err_full[i, :4])
        s02_val, s02_err = convert_to_str(best_Fit_r[i][0], prec[0]), convert_to_str(err_full[i, 0], prec[0])
        e0_val, e0_err = convert_to_str(best_Fit_r[i][1], prec[1]), convert_to_str(err_full[i, 1], prec[1])
        sig_val, sig_err = np.round(best_Fit_r[i][2], 4), np.round(err_full[i, 2], 4)
        r_val = np.round(og_arr_r[i], 4)
        fr_val, fr_err = np.round(best_Fit_r[i][3], 3), temp_round_deltaR(err_full[i, 3])
        legs = int(best_Fit_r[i][5])

        latex_table_str += (
            f"                        {p} & {int(best_Fit_r[i][4])} & {s02_val}$\\pm${s02_err} & "
            f"{e0_val}$\\pm${e0_err} & {sig_val}$\\pm${sig_err} & {r_val} & {fr_val}$\\pm${fr_err} & {legs} & {label_arr[i]}\\\\"
            "\n"
        )

    latex_table_str += """                        \\hline
                        \\end{tabular}
                        \\label{Label}
\\end{table}"""

    return nleg_arr, label_arr, latex_table_str