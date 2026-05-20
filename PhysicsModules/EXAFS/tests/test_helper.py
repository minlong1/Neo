import unittest

from PhysicsModules.EXAFS.exafs_neo.helper import str_to_bool


class TestValidator(unittest.TestCase):
    def test_str2bool_true(self):
        self.assertEqual(True,str_to_bool("True"))

    def test_str2bool_false(self):
        self.assertEqual(False,str_to_bool('False'))

if __name__ == '__main__':
    unittest.main()
