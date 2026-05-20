import logging
import unittest

from PhysicsModules.EXAFS.exafs_neo.utils import NeoLogger


class TestUtils(unittest.TestCase):

    def test_logger_level(self):
        exafs_logger = NeoLogger()
        exafs_logger.initialize_logging()
        self.assertEqual(exafs_logger.logging_level, logging.DEBUG)

    def test_checkoutput_file(self):
        pass
