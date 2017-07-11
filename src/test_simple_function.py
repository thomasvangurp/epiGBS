import unittest
import simple_function

class IntegerArithmeticTestCase(unittest.TestCase):
    def testAdd(self): ## test method names begin 'test*'
        self.assertEqual(simple_function.add_a_b(2,3),6)


if __name__ == '__main__':
    unittest.main()