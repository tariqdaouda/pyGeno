import unittest
import sys
import os

_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _root not in sys.path:
    sys.path.insert(0, _root)

import numpy as np
from tools.Stats import kullback_leibler, squaredError_log10, fisherExactTest


class TestKullbackLeibler(unittest.TestCase):

    def test_identical_distributions(self):
        p = [0.25, 0.25, 0.25, 0.25]
        q = [0.25, 0.25, 0.25, 0.25]
        result = kullback_leibler(p, q)
        self.assertAlmostEqual(result, 0.0, places=5)

    def test_different_distributions(self):
        p = [0.5, 0.5]
        q = [0.25, 0.75]
        result = kullback_leibler(p, q)
        self.assertGreater(result, 0)

    def test_kl_nonnegative(self):
        p = [0.1, 0.9]
        q = [0.5, 0.5]
        result = kullback_leibler(p, q)
        self.assertGreaterEqual(result, 0)

    def test_mismatched_shapes_raises(self):
        p = [0.5, 0.5]
        q = [0.25, 0.25, 0.5]
        with self.assertRaises(ValueError):
            kullback_leibler(p, q)

    def test_zero_in_p_handled(self):
        p = [0.0, 1.0]
        q = [0.5, 0.5]
        result = kullback_leibler(p, q)
        # log(1.0/0.5) * 1.0 = log(2)
        self.assertAlmostEqual(result, np.log(2), places=5)

    def test_returns_float(self):
        p = [0.5, 0.5]
        q = [0.5, 0.5]
        result = kullback_leibler(p, q)
        self.assertIsInstance(float(result), float)


class TestSquaredErrorLog10(unittest.TestCase):

    def test_identical_arrays(self):
        p = [1.0, 2.0, 3.0]
        q = [1.0, 2.0, 3.0]
        result = squaredError_log10(p, q)
        # sum of (p-q)^2 = 0, log10(0) = -inf
        self.assertEqual(result, float('-inf') - np.log(3))

    def test_different_arrays(self):
        p = [1.0, 2.0]
        q = [2.0, 3.0]
        result = squaredError_log10(p, q)
        # (1-2)^2 + (2-3)^2 = 2, log10(2) - ln(2)
        expected = np.log10(2) - np.log(2)
        self.assertAlmostEqual(float(result), float(expected), places=5)

    def test_mismatched_shapes_raises(self):
        p = [1.0, 2.0]
        q = [1.0, 2.0, 3.0]
        with self.assertRaises(ValueError):
            squaredError_log10(p, q)


class TestFisherExact(unittest.TestCase):

    def test_not_implemented(self):
        with self.assertRaises(NotImplementedError):
            fisherExactTest([[1, 2], [3, 4]])


class TestNumpyDtypeCompatibility(unittest.TestCase):
    """Verify the np.float -> float fix works."""

    def test_kl_uses_float_dtype(self):
        p = [0.5, 0.5]
        q = [0.5, 0.5]
        # This would fail with np.float on NumPy >= 1.24
        result = kullback_leibler(p, q)
        self.assertIsNotNone(result)

    def test_squared_error_uses_float_dtype(self):
        p = [1.0, 2.0]
        q = [1.0, 2.0]
        result = squaredError_log10(p, q)
        self.assertIsNotNone(result)


if __name__ == "__main__":
    unittest.main()
