#! /usr/bin/env python

from __future__ import print_function
from . import get_eb_score
import sys, os, subprocess, math, re, multiprocessing 
import pysam, numpy

region_exp = re.compile('^([^ \t\n\r\f\v,]+):(\d+)\-(\d+)')

def EBFilter_worker_vcf(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode, control_database):

    import vcfpy
    from . import process_vcf

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))

    ##########
    # generate pileup files
    process_vcf.vcf2pileup(targetMutationFile, outputPath + '.target.pileup', targetBamPath, mapping_qual_thres, base_qual_thres, filter_flags, False, is_loption, region)

    ##########
    # load pileup files
    pos2pileup_target = {}

    hIN = open(outputPath + '.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    ##########
    # get restricted region if not None
    if is_loption == True and region != "":
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))

    ##########
    # control database
    ctrl_db_tabix = pysam.TabixFile(control_database)

    vcf_reader = vcfpy.Reader.from_path(targetMutationFile)

    vcf_reader.header.add_info_line(vcfpy.OrderedDict(
        [('ID','EB'), ('Number','1'), ('Type','Float'), ('Description','EBCall Score')]))
    vcf_writer = vcfpy.Writer.from_path(outputPath, vcf_reader.header)

    for vcf_record in vcf_reader:
        current_pos = str(vcf_record.CHROM) + '\t' + str(vcf_record.POS) 

        if is_loption == True and region != "":
            if reg_chr != vcf_record.CHROM: continue
            if int(vcf_record.POS) < reg_start or int(vcf_record.POS) > reg_end: continue

        F_target = pos2pileup_target[current_pos].split('\t') if current_pos in pos2pileup_target else []

        current_ref = str(vcf_record.REF)
        current_alt = str(vcf_record.ALT[0].value)
        var = ""
        if len(current_ref) == 1 and len(current_alt) == 1:
            var = current_alt
        else:
            if len(current_ref) == 1:
                var = "+" + current_alt[1:]
            elif len(current_alt) == 1:
                var = "-" + current_ref[1:]

        EB_score = "." # if the variant is complex, we ignore that
        if not var == "":

            records = None
            param_region = f"{vcf_record.CHROM}:{vcf_record.POS}-{vcf_record.POS}"
            try:
                records = ctrl_db_tabix.fetch(region=param_region)
            except Exception as inst:
                print(inst.args, file = sys.stderr)

            alpha_p, beta_p, alpha_n, beta_n = None, None, None, None
            if records != None:
                for record_line in records:
                    F = record_line.split("\t")
                    if F[0] == vcf_record.CHROM and F[1] == str(vcf_record.POS):
                        if ((F[4] == var)
                        or (var.startswith("+") and F[4] == "<INS>")
                        or (var.startswith("-") and F[4] == "<DEL>")):
                            # debug
                            # print(f"Using Database! {vcf_record.CHROM},{vcf_record.POS},{vcf_record.REF},{vcf_record.ALT[0].value}")
                            for item in F[7].split(";"):
                                key,val = item.split("=")
                                if key == "AP": alpha_p = numpy.float64(val)
                                if key == "BP": beta_p = numpy.float64(val)
                                if key == "AN": alpha_n = numpy.float64(val)
                                if key == "BN": beta_n = numpy.float64(val)

            if alpha_p == None or beta_p == None or alpha_n == None or beta_n == None:
                # debug
                # print(f"NOT Using Database! {vcf_record.CHROM},{vcf_record.POS},{vcf_record.REF},{vcf_record.ALT[0].value}")
                F_control = process_vcf.pileup1line(controlBamPathList, mapping_qual_thres, base_qual_thres, filter_flags, param_region)
                alpha_p, beta_p, alpha_n, beta_n = get_eb_score.get_beta_binomial_alpha_and_beta(var, F_control, base_qual_thres, controlFileNum)

            EB_score = get_eb_score.get_eb_score(var, F_target, base_qual_thres, controlFileNum, alpha_p, beta_p, alpha_n, beta_n)

        # add the score and write the vcf record
        vcf_record.INFO['EB'] = EB_score
        vcf_writer.write_record(vcf_record)

    vcf_writer.close()


    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", outputPath + '.target.pileup'])


def EBFilter_worker_anno(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode, control_database):

    from . import process_anno

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))

    ##########
    # generate pileup files
    process_anno.anno2pileup(targetMutationFile, outputPath + '.target.pileup', targetBamPath, mapping_qual_thres, base_qual_thres, filter_flags, False, is_loption, region)

    ##########
    # load pileup files
    pos2pileup_target = {}

    hIN = open(outputPath + '.target.pileup')
    for line in hIN:
        F = line.rstrip('\n').split('\t')
        pos2pileup_target[F[0] + '\t' + F[1]] = '\t'.join(F[3:])
    hIN.close()

    ##########
    # get restricted region if not None
    if is_loption == True and region != "":
        region_match = region_exp.match(region)
        reg_chr = region_match.group(1)
        reg_start = int(region_match.group(2))
        reg_end = int(region_match.group(3))

    ctrl_db_tabix = pysam.TabixFile(control_database)
    hIN = open(targetMutationFile, 'r')
    hOUT = open(outputPath, 'w')

    for line in hIN:

        inF = line.rstrip('\n').split('\t')
        chrom, pos, pos2, ref, alt = inF[0], inF[1], inF[2], inF[3], inF[4]
        if alt == "-": pos = str(int(pos) - 1)

        if is_loption == True and region != "":
            if reg_chr != chrom: continue
            if int(pos) < reg_start or int(pos) > reg_end: continue

        F_target = pos2pileup_target[chrom + '\t' + pos].split('\t') if chrom + '\t' + pos in pos2pileup_target else []

        var = ""
        if ref != "-" and alt != "-":
            var = alt
        else:
            if ref == "-":
                var = "+" + alt
            elif alt == "-":
                var = "-" + ref

        EB_score = "." # if the variant is complex, we ignore that
        if not var == "":

            records = None
            param_region = f"{chrom}:{pos}-{pos}"
            try:
                records = ctrl_db_tabix.fetch(region=param_region)
            except Exception as inst:
                print(inst.args, file = sys.stderr)

            alpha_p, beta_p, alpha_n, beta_n = None, None, None, None
            
            if records != None:
                for record_line in records:
                    F = record_line.split("\t")
                    if F[0] == chrom and F[1] == pos:
                        if ((F[4] == var)
                        or (var.startswith("+") and F[4] == "<INS>")
                        or (var.startswith("-") and F[4] == "<DEL>")):
                            # debug
                            # print(f"Using Database! {chrom},{pos},{ref},{alt}")
                            for item in F[7].split(";"):
                                key,val = item.split("=")
                                if key == "AP": alpha_p = numpy.float64(val)
                                if key == "BP": beta_p = numpy.float64(val)
                                if key == "AN": alpha_n = numpy.float64(val)
                                if key == "BN": beta_n = numpy.float64(val)

            if alpha_p == None or beta_p == None or alpha_n == None or beta_n == None:
                # debug
                # print(f"NOT Using Database! {chrom},{pos},{ref},{alt}")
                F_control = process_anno.pileup1line(controlBamPathList, mapping_qual_thres, base_qual_thres, filter_flags, param_region)
                alpha_p, beta_p, alpha_n, beta_n = get_eb_score.get_beta_binomial_alpha_and_beta(var, F_control, base_qual_thres, controlFileNum)

            EB_score = get_eb_score.get_eb_score(var, F_target, base_qual_thres, controlFileNum, alpha_p, beta_p, alpha_n, beta_n)

        # add the score and write the vcf record
        print('\t'.join(inF + [str(EB_score)]), file=hOUT)

    hIN.close()
    hOUT.close()


    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", outputPath + '.target.pileup'])


def ebfilter_main(args):

    # should add validity check for arguments
    targetMutationFile = args.targetMutationFile
    targetBamPath = args.targetBamPath
    controlBamPathList = args.controlBamPathList
    outputPath = args.outputPath

    mapping_qual_thres = args.q
    base_qual_thres = args.Q
    filter_flags = args.ff
    thread_num = args.t
    is_anno = True if args.f == 'anno' else False
    is_loption = args.loption
    region = args.region
    debug_mode = args.debug
    control_panel_database = args.panel

    # region format check
    if region != "":
        region_match = region_exp.match(region)
        if region_match is None:
            print("Wrong format for --region ({chr}:{start}-{end}): " + region, file=sys.stderr)
            sys.exit(1)

    # file existence check
    if not os.path.exists(targetMutationFile):
        print("No target mutation file: " + targetMutationFile, file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(targetBamPath):
        print("No target bam file: " + targetBamPath, file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(targetBamPath + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", targetBamPath)):
        print("No index for target bam file: " + targetBamPath, file=sys.stderr)
        sys.exit(1)


    if not os.path.exists(controlBamPathList):
        print("No control list file: " + controlBamPathList, file=sys.stderr)
        sys.exit(1)

    with open(controlBamPathList) as hIN:
        for in_file in hIN:
            in_file = in_file.rstrip()
            if not os.path.exists(in_file):
                print("No control bam file: " + in_file, file=sys.stderr) 
                sys.exit(1)

            if not os.path.exists(in_file + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", in_file)):
                print("No index control bam file: " + in_file, file=sys.stderr) 
                sys.exit(1)

    if thread_num == 1:
        # non multi-threading mode
        if is_anno == True:
            EBFilter_worker_anno(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode, control_panel_database)
        else: 
            EBFilter_worker_vcf(targetMutationFile, targetBamPath, controlBamPathList, outputPath, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode, control_panel_database)
    else:
        # multi-threading mode
        ##########

        if is_anno == True:
            from . import process_anno
            # partition anno files
            process_anno.partition_anno(targetMutationFile, outputPath + ".tmp.input.anno.", thread_num)

            jobs = []
            for i in range(thread_num):
                process = multiprocessing.Process(target = EBFilter_worker_anno, args = \
                    (outputPath + ".tmp.input.anno." + str(i), targetBamPath, controlBamPathList, outputPath + "." + str(i), mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode, control_panel_database))
                jobs.append(process)
                process.start()
        
            # wait all the jobs to be done
            for i in range(thread_num):
                jobs[i].join()
        
            # merge the individual results
            process_anno.merge_anno(outputPath + ".", outputPath, thread_num)
        
            # delete intermediate files
            if debug_mode == False:
                for i in range(thread_num):
                    subprocess.check_call(["rm", outputPath + ".tmp.input.anno." + str(i)])
                    subprocess.check_call(["rm", outputPath + "." + str(i)])

        else:
            from . import process_vcf
            # partition vcf files
            process_vcf.partition_vcf(targetMutationFile, outputPath + ".tmp.input.vcf.", thread_num)

            jobs = []
            for i in range(thread_num):
                process = multiprocessing.Process(target = EBFilter_worker_vcf, args = \
                    (outputPath + ".tmp.input.vcf." + str(i), targetBamPath, controlBamPathList, outputPath + "." + str(i), mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode, control_panel_database))
                jobs.append(process)
                process.start()

            # wait all the jobs to be done
            for i in range(thread_num):
                jobs[i].join()

            # merge the individual results
            process_vcf.merge_vcf(outputPath + ".", outputPath, thread_num)

            # delete intermediate files
            if debug_mode == False:
                for i in range(thread_num):
                    subprocess.check_call(["rm", outputPath + ".tmp.input.vcf." + str(i)])
                    subprocess.check_call(["rm", outputPath + "." + str(i)])


def create_control_panel_database_worker(targetMutationFile, controlBamPathList, outputPathPrefix, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode):

    import vcfpy
    from . import process_vcf

    cur_chrom = None
    cur_pos = None

    controlFileNum = sum(1 for line in open(controlBamPathList, 'r'))

    with open(targetMutationFile, "r") as hin, open(outputPathPrefix + '.tmp.vcf', "w") as hout:
        for line in hin:
            if line.startswith("#"):
                hout.write(line)
                continue

            F = line.split("\t")
            if cur_chrom is None or F[0] != cur_chrom or F[1] != cur_pos:
                hout.write(line)

            cur_chrom = F[0]
            cur_pos = F[1]


    ##########
    # generate pileup files
    process_vcf.vcf2pileup(outputPathPrefix +'.tmp.vcf', outputPathPrefix + '.control.pileup', controlBamPathList, mapping_qual_thres, base_qual_thres, filter_flags, True, is_loption, region)

    vcf_reader = vcfpy.Reader.from_path(targetMutationFile)

    with open(outputPathPrefix + ".vcf", "w") as hout, open(outputPathPrefix + '.control.pileup', 'r') as hin:
        hout.write('##fileformat=VCFv4.2\n')
        hout.write('##reference=GRCh38\n')
        hout.write('##INFO=<ID=AP,Number=1,Type=Float,Description="Alpha which estimates the beta-binomial parameters for positive strands">\n')
        hout.write('##INFO=<ID=BP,Number=1,Type=Float,Description="Beta which estimates the beta-binomial parameters for positive strands">\n')
        hout.write('##INFO=<ID=AN,Number=1,Type=Float,Description="Alpha which estimates the beta-binomial parameters for negative strands">\n')
        hout.write('##INFO=<ID=BN,Number=1,Type=Float,Description="Beta which estimates the beta-binomial parameters for negative strands">\n')
        hout.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

        pileup = hin.readline()
        l_pileup = pileup.strip("\n").split("\t") if pileup else None

        # count = 0
        for vcf_record in vcf_reader:

            while l_pileup and (l_pileup[0] != str(vcf_record.CHROM) or (l_pileup[0] == str(vcf_record.CHROM) and int(l_pileup[1]) < int(vcf_record.POS))):
                pileup = hin.readline()
                l_pileup = pileup.strip("\n").split("\t") if pileup else None
            F_control = []
            if l_pileup and (l_pileup[0] == str(vcf_record.CHROM) and int(l_pileup[1]) == int(vcf_record.POS)):
                F_control = l_pileup[3:]

            current_ref = str(vcf_record.REF)
            current_alt = str(vcf_record.ALT[0].value)
            var = ""
            if len(current_ref) == 1 and len(current_alt) == 1:
                var = current_alt
            elif current_alt == "INS":
                var = "+N"
                current_alt = "<INS>"
            elif current_alt == "DEL":
                var = "-N"
                current_alt = "<DEL>"

            alpha_p, beta_p, alpha_n, beta_n = get_eb_score.get_beta_binomial_alpha_and_beta(var, F_control, base_qual_thres, controlFileNum)

            # add the score and write the vcf record
            hout.write(f"{str(vcf_record.CHROM)}\t{vcf_record.POS}\t.\t{current_ref}\t{current_alt}\t.\t.\tAP={round(alpha_p, 4)};BP={round(beta_p, 4)};AN={round(alpha_n, 4)};BN={round(beta_n, 4)}\n")

            # count += 1
            # if count % 10000 == 0: print(f"vcf line: {count}")

    # delete intermediate files
    if debug_mode == False:
        subprocess.check_call(["rm", outputPathPrefix + '.control.pileup'])
        subprocess.check_call(["rm", outputPathPrefix + '.tmp.vcf'])

    subprocess.check_call(["bgzip", "-f", outputPathPrefix + ".vcf"])
    subprocess.check_call(["tabix", "-p", "vcf", outputPathPrefix + ".vcf.gz"])


def create_control_panel_database(args):
    targetMutationFile = args.targetMutationFile
    controlBamPathList = args.controlBamPathList
    outputPathPrefix = args.outputPathPrefix
    is_loption = args.loption
    region = args.region
    mapping_qual_thres = args.q
    base_qual_thres = args.Q
    filter_flags = args.ff
    debug_mode = args.debug

    # region format check
    if region != "":
        region_match = region_exp.match(region)
        if region_match is None:
            print("Wrong format for --region ({chr}:{start}-{end}): " + region, file=sys.stderr)
            sys.exit(1)

    # file existence check
    if not os.path.exists(targetMutationFile):
        print("No target mutation file: " + targetMutationFile, file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(controlBamPathList):
        print("No control list file: " + controlBamPathList, file=sys.stderr)
        sys.exit(1)

    with open(controlBamPathList) as hIN:
        for in_file in hIN:
            in_file = in_file.rstrip()
            if not os.path.exists(in_file):
                print("No control bam file: " + in_file, file=sys.stderr)
                sys.exit(1)

            if not os.path.exists(in_file + ".bai") and not os.path.exists(re.sub(r'bam$', "bai", in_file)):
                print("No index control bam file: " + in_file, file=sys.stderr)
                sys.exit(1)

    create_control_panel_database_worker(targetMutationFile, controlBamPathList, outputPathPrefix, mapping_qual_thres, base_qual_thres, filter_flags, is_loption, region, debug_mode)


