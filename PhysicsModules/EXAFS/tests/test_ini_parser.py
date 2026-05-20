import unittest

from PhysicsModules.EXAFS.exafs_neo.ini_parser import optional_var, split_path_arr


class TestOptionalVar(unittest.TestCase):
    def test_optional_var(self):
        """
        Testing if it is the correct value and type
        :return:
        """
        input_dict = {}
        optional_var(input_dict, 'file', 0, int)
        self.assertEqual(input_dict['file'], 0)
        self.assertIsInstance(input_dict['file'], int)

    def test_optional_bool(self):
        """
        Testing for boolean
        :return:
        """
        input_dict = {}
        optional_var(input_dict, 'file', True, bool)
        self.assertEqual(input_dict['file'], True)
        self.assertIsInstance(input_dict['file'], bool)


class TestSplitPathArr(unittest.TestCase):
    def test_split_path_arr_single(self):
        """
        Test Single components array
        :return:
        """
        input_arr = '1,2,3'
        result_str = split_path_arr(input_arr, 1)
        self.assertEqual(['1', '2', '3'], result_str)

    def test_split_path_arr_single_bracket(self):
        """
        Test Single components array
        :return:
        """
        input_arr = '[1,2,3]'
        result_str = split_path_arr(input_arr, 1)
        self.assertEqual([['1', '2', '3']], result_str)

    def test_split_path_arr_two_comp(self):
        """
        Test Single components array
        :return:
        """
        input_arr = '[1,2,3],[1,2]'
        result_str = split_path_arr(input_arr, 1)
        self.assertEqual([['1', '2', '3'], ['1', '2']], result_str)


if __name__ == '__main__':
    unittest.main()
