#! /usr/bin/env python

from __future__ import print_function 
import unittest
import os, tempfile, shutil, filecmp
import subprocess

class TestEBFilter(unittest.TestCase):

    def test1(self):
        subprocess.check_call(["EBFilter", "--version"])

    def test2(self):
        subprocess.check_call(["samtools", "--version"])

    def test3(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_anno = cur_dir + "/../testdata/input.anno"
        in_bam = cur_dir + "/../testdata/tumor.bam"
        in_panel = cur_dir + "/../testdata/control_panel1.vcf.gz"
        output_file = tmp_dir+"/output.anno"
        control_panel = tmp_dir + "/list_normal_sample.txt"
        with open (control_panel, "w") as hout:
            print(cur_dir + "/../testdata/normalreference1.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference2.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference3.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference4.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference5.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference6.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference7.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference8.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference9.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference10.bam", file=hout)

        subprocess.check_call(["EBFilter", "main", "-f", "anno", in_anno, in_bam, output_file, "--panel", in_panel, "-c", control_panel])

        answer_file = cur_dir + "/../testdata/output.golden.anno"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test4(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_vcf = cur_dir + "/../testdata/input.vcf.gz"
        in_bam = cur_dir + "/../testdata/tumor.bam"
        in_panel = cur_dir + "/../testdata/control_panel2.vcf.gz"
        output_file = tmp_dir+"/output.vcf"
        control_panel = tmp_dir + "/list_normal_sample.txt"
        with open (control_panel, "w") as hout:
            print(cur_dir + "/../testdata/normalreference1.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference2.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference3.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference4.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference5.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference6.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference7.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference8.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference9.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference10.bam", file=hout)

        subprocess.check_call(["EBFilter", "main", "-f", "vcf", in_vcf, in_bam, output_file, "--panel", in_panel, "-c", control_panel])

        answer_file = cur_dir + "/../testdata/output.golden.vcf"
        self.assertTrue(filecmp.cmp(output_file, answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

    def test5(self):
        cur_dir = os.path.dirname(os.path.abspath(__file__))
        tmp_dir = tempfile.mkdtemp()
        in_vcf = cur_dir + "/../testdata/input.vcf"
        output_prefix = tmp_dir+"/control_panel"
        control_panel = tmp_dir + "/list_normal_sample.txt"
        with open (control_panel, "w") as hout:
            print(cur_dir + "/../testdata/normalreference1.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference2.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference3.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference4.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference5.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference6.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference7.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference8.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference9.bam", file=hout)
            print(cur_dir + "/../testdata/normalreference10.bam", file=hout)

        subprocess.check_call(["EBFilter", "panel", in_vcf, control_panel, output_prefix])

        answer_file = cur_dir + "/../testdata/control_panel_golden.vcf.gz"
        self.assertTrue(filecmp.cmp(output_prefix+".vcf.gz", answer_file, shallow=False))
        shutil.rmtree(tmp_dir)

