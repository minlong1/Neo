import numpy as np


class xps_data:
    """Spectrum loader, unified from the CLI and GUI forks (Phase 3).

    Divergence decisions (docs/fork-divergence.md):
    - Offsets are applied ONCE at load. The CLI fork re-applied them on
      every get_x/get_y call (double-applying y_offset per fit setup) and
      the GUI fork applied them only inside noise_check/range_check; both
      only ever ran with offset 0.0 in golden-covered paths, where all
      three behaviors coincide.
    - get_x takes the GUI's (data_KE, data_XES) flags with CLI-compatible
      defaults; get_y and the scale logic were already identical.
    - noise_check/range_check come from the GUI fork (minus their offset
      side effects, now handled at load).
    """

    def __init__(self, fileName, skipLn, x_offset, y_offset, data_KE=False):
        self.fileName = fileName
        self.skipLn = skipLn  # int represnting files to skip
        self.x_offset = x_offset
        self.y_offset = y_offset
        self.read_file_txt()

    def read_file_txt(self, fileName=None):
        self.x, self.y = np.loadtxt(self.fileName, skiprows=self.skipLn, unpack=True)
        self.numP = len(self.x)
        for i in range(len(self.x)):
            self.x[i] = self.x[i] + self.x_offset
        for i in range(len(self.y)):
            self.y[i] = self.y[i] + self.y_offset

    def noise_check(self):
        # check to see if data is too noisy
        # finds the average on each side of data and sees how close it is to maximum y value
        N = 30
        y_left = self.y[:N]
        y_right = self.y[-N:]
        y_left_avg = sum(y_left) / N
        y_right_avg = sum(y_right) / N
        y_max = np.max(self.y)
        # if +/- 10% of average of sides --> Too low?
        if y_left_avg <= y_max <= y_left_avg + 0.1 * y_left_avg:
            print("Warning: Data is too noisy to be fit. Results will be moot")
            print("To continue analysis select 'Enable expert mode' in the Expert tab")
            return False
        elif y_right_avg <= y_max <= y_right_avg + 0.1 * y_right_avg:
            print("Warning: Data is too noisy to be fit. Results will be moot")
            print("To continue analysis select 'Enable expert mode' in the Expert tab")
            return False
        else:
            return True

    def range_check(self):
        # check if data has wide enough range for fitting
        numP = len(self.y)
        y_max = np.max(self.y)

        for i in range(numP):
            if self.y[i] == y_max:
                x_max = self.x[i]
        x_left = self.x[0]
        x_right = self.x[-1]
        if numP < 100:  # What is an appropriate number of data points?
            print("Warning: Data range is too small. Results will be moot")
            print("To continue analysis select 'Enable expert mode' in the Expert tab")

        # Is +/- 5 eV good?
        if x_left > x_right:
            if x_left - 5 <= x_max <= x_left:
                print("Warning: Data range is too small. Results will be moot")
                print(
                    "To continue analysis select 'Enable expert mode' in the Expert tab"
                )
                return False
            elif x_right <= x_max <= x_right + 5:
                print("Warning: Data range is too small. Results will be moot")
                print(
                    "To continue analysis select 'Enable expert mode' in the Expert tab"
                )
                return False
            else:
                return True
        else:  # KE
            if x_left <= x_max <= x_left + 5:
                print("Warning: Data range is too small. Results will be moot")
                print(
                    "To continue analysis select 'Enable expert mode' in the Expert tab"
                )
                return False
            elif x_right - 5 <= x_max <= x_right:
                print("Warning: Data range is too small. Results will be moot")
                print(
                    "To continue analysis select 'Enable expert mode' in the Expert tab"
                )
                return False
            else:
                return True

    def get_x(self, data_KE=False, data_XES=False):
        if data_KE == True:
            self.x = self.x[::-1]
        elif data_XES == True:
            self.x = self.x[::-1]
        return self.x

    def get_y(self, scale_var):
        y_val = []
        first = self.y[0]
        if first < 1:
            scale_val = round(1 / first) * 100
        else:
            scale_val = first

        if scale_var == True:
            y_val = self.y / first  # Dividing every element by the first value
            y_val = self.y * scale_val  # Multiply by 1000 to scale it
        else:
            y_val = self.y
        return y_val
