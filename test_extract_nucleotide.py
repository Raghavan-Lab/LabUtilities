# coding=utf-8
import unittest
import extract_nucleotide


class TestExtractNucleotide(unittest.TestCase):
    """ Simple unit test to verify the output of the extract nucleotide script.  Known-good sequences from
    NCBI are compared against the script results """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def test_NP_820224_nucleotides(self):
        # build expected results
        expected_result_file = '/Users/tsmit2/Dropbox/PycharmProjects/RaghavanLab/GitHub/LabUtilities/test_resources/NP_820224_expected_nucleotide_sequence.fa'
        expected_result = open(expected_result_file).read()

        test_list = 'test_resources/coxiella_test_NP_820224.txt'
        genbank_file = 'test_resources/Coxiella.gb'

        for result in extract_nucleotide.main(genbank_file=genbank_file,
                                              target_list=test_list):
            print('Expect:')
            print(expected_result)
            print('---------------')
            print('Produced:')
            print(result)
            self.assertEqual(expected_result.rstrip(), result.rstrip(), msg='Output does not match expected result')


if __name__ == '__main__':
    unittest.main()