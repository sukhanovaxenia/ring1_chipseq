import importlib.util
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('pipeline',ROOT/'scripts/run.py')
pipeline=importlib.util.module_from_spec(spec);spec.loader.exec_module(pipeline)

class BedTests(unittest.TestCase):
    def read(self,text):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'x.bed';p.write_text(text)
            return pipeline.read_bed(p)
    def test_real_bdgdiff_bed5_and_metadata(self):
        x=self.read('track name=test\n# comment\nchr1\t16535157\t16535414\tring/peaks/diff/RING1_BF_unique.bed1\t9.40754\n')
        self.assertEqual(x,[('chr1',16535157,16535414)])
    def test_invalid_coordinates(self):
        for text in ['chr1 -1 8','chr1 9 9','chr1 10 9','chr1 x 8','chr1 4']:
            with self.subTest(text=text),self.assertRaises(ValueError):self.read(text)
    def test_narrowpeak_first_three_fields(self):
        self.assertEqual(self.read('chrX 0 10 p 12 . 5 6 7 2'),[('chrX',0,10)])
    def test_empty_file(self):self.assertEqual(self.read(''),[])
    def test_union_dedup_touching_nested_and_chromosomes(self):
        rows=[('chr2',0,10),('chr1',3,5),('chr1',0,10),('chr1',10,15),('chr1',0,10),('chr1',20,21)]
        self.assertEqual(pipeline.merge_intervals(rows),[('chr1',0,15),('chr1',20,21),('chr2',0,10)])
    def test_union_preserves_gaps(self):
        self.assertEqual(pipeline.merge_intervals([('chr1',0,5),('chr1',6,10)]),[('chr1',0,5),('chr1',6,10)])
    def test_bed3_roundtrip(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'roundtrip.bed';rows=[('chr1',0,10),('chrX',12,14)]
            pipeline.write_bed(p,rows)
            self.assertEqual(pipeline.read_bed(p),rows)
            self.assertTrue(all(len(line.split('\t'))==3 for line in p.read_text().splitlines()))

if __name__=='__main__':unittest.main()
