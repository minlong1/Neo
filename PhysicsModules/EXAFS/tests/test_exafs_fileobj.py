import unittest

from PhysicsModules.EXAFS.exafs_neo.exafsfileobj import ExafsFileobj


class MyTestCase(unittest.TestCase):

    def test_neoPars_single(self):
        exafs_File = ExafsFileobj()

        data_dict = {
            "data_file": "path_files/Pu_C/pu3in_t030_sac.chi",
            "feff_file": "path_files/Pu_C/feff/feff",
            "output_file": "result/test/PuC.csv",
            "sabcor_toggle": False
        }

        exafs_File.initialize(data_dict)

        self.assertEqual(exafs_File.csv_series, False)
        self.assertEqual(exafs_File.data_file,'path_files/Pu_C/pu3in_t030_sac.chi')
        self.assertEqual(exafs_File.end,'.dat')
        self.assertEqual(exafs_File.nComp,1)
        self.assertEqual(exafs_File.sabcorFile, False)
        self.assertEqual(exafs_File.sabcor_toggle, False)

    def test_neoPars_multiple(self):
        exafs_File = ExafsFileobj()

        data_dict = {
            "data_file": "path_files/Pu_C/pu3in_t030_sac.chi",
            "feff_file": ["path_files/Pu_C/feff/feff", "path_files/Pu_C/feff/feff2"],
            "output_file": "result/test/PuC.csv",
            "nComp": 2,
            "sabcor_toggle": False
        }

        exafs_File.initialize(data_dict)

        self.assertEqual(exafs_File.feff_file[0],'path_files/Pu_C/feff/feff')
        self.assertEqual(exafs_File.feff_file[1],'path_files/Pu_C/feff/feff2')
        self.assertEqual(exafs_File.nComp,2)


    def test_neoPars_sabcor(self):
        pass


    def test_neoPars_check_outputFile(self):
        pass
# add assertion here


if __name__ == '__main__':
    unittest.main()
