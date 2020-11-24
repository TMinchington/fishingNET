"""
fishingNet is a suite of tools for extracting data from the network
"""


def get_regulator_targets(network_file, gene, filterX, outpath):

    """
    Outputs all of the targets of the genes in "gene" in the "network_file"
    If a filterX is not false targets will be filtered by the target type provided
    """

    if not filterX:
        outfile = open(os.path.join(outpath, f'{"-".join(gene)}_targets.tsv'), 'w')

    else:
        outfile = open(os.path.join(outpath, f'{"-".join(gene)}_{filterX}_targets.tsv'), 'w')

    with open(network_file) as openNet:
        first_line = True
        for line in openNet:
            split_line = line.strip().split('\t')
            if first_line:
                outfile.write(line)
                genNameDex = split_line.index('regulator_name')
                targetTypeDex = split_line.index('target_type')
                first_line = False
                continue
            
            geneName = split_line[genNameDex]
            targetType = split_line[targetTypeDex]

            for x in gene:
                if x == geneName:
                    if filterX:
                        if targetType == filterX:
                            outfile.write(line)
                            break
                        continue
                    outfile.write(line)
                    break
    
    outfile.close()

    


def get_target_regulators(network_file, gene, filterX, outpath):

    """
    Outputs all of the regulators of the genes in "gene" in the "network_file"
    If a filterX is not false, regualtors will be filtered by the regulator type provided
    """

    if not filterX:
        outfile = open(os.path.join(outpath, f'{"-".join(gene)}_regulators.tsv'), 'w')

    else:
        outfile = open(os.path.join(outpath, f'{"-".join(gene)}_{filterX}_regulators.tsv'), 'w')

    with open(network_file) as openNet:
        first_line = True
        for line in openNet:
            split_line = line.strip().split('\t')
            if first_line:
                outfile.write(line)
                genNameDex = split_line.index('target_name')
                regTypeDex = split_line.index('reg_type')
                first_line = False
                continue
            
            geneName = split_line[genNameDex]
            regType = split_line[regTypeDex]

            for x in gene:
                if x == geneName:
                    if filterX:
                        if regType == filterX:
                            outfile.write(line)
                            break
                        continue
                    outfile.write(line)
                    break
    
    outfile.close()


def surroundingNetwork(network_file, gene, outpath):

    """
    Finds all the inputs and outputs of one gene
    """
    outfile = open(os.path.join(outpath, f'{"-".join(gene)}_regulators_and_target.tsv'), 'w')

    with open(network_file) as openNet:
        first_line = True
        for line in openNet:
            split_line = line.strip().split('\t')
            if first_line:
                outfile.write(line)
                regNameDex = split_line.index('target_name')
                tarNameDex = split_line.index('regulator_name')
                first_line = False
                continue
            
            regName = split_line[regNameDex]
            tarName = split_line[tarNameDex]

            for x in gene:
                if x == regName or x == tarName:
                    outfile.write(line)
    
    outfile.close()    


def checkCrCt(network_file, cr, ct, outpath):

    """
    Checks to see if two partial genes interact
    """

    cr = cr.lower()
    ct = ct.lower()

    outfile = open(os.path.join(outpath, f'{cr}_{ct}_interaction.tsv'), 'w')
    counter = 0
    count_dic = {}
    with open(network_file) as openNet:
        first_line = True
        for line in openNet:
            split_line = line.strip().split('\t')
            if first_line:
                outfile.write(line)
                tarNameDex = split_line.index('target_name')
                regNameDex = split_line.index('regulator_name')
                first_line = False
                continue
            
            regName = split_line[regNameDex].lower()
            tarName = split_line[tarNameDex].lower()

            if cr in regName and ct in tarName:
                counter += 1
                try:
                    count_dic[regName] += 1

                except KeyError:
                    count_dic[regName] = 1

                try:
                    count_dic[tarName] += 1
                except KeyError:
                    count_dic[tarName] = 1

                outfile.write(line)
    
    print('\n\n-----------------------------------------------\n\n')

    if counter == 0:
        print(f"No matches found for {cr} regulating {ct}")

    else:
        print(f"Matches found {cr}:\n")
        for x in count_dic:
            if cr in x:
                print('\t\t', x, count_dic[x])

        print(f"\nMatches for {ct}:\n")
        for x in count_dic:
            if ct in x:
                print('\t\t', x, count_dic[x])

    print('\n\n-----------------------------------------------\n\n')
    outfile.close()    

def random_snarky_comment(checkMe):
    from random import randint
    ran_ls = ["HES genes again how original", "Don't we search for HES every week?", "HES! AGAIN! Come on....", "HES? Yep sure why not, I can search for anything though you know?",
     "HES gon' give it to ya, hes gon' give it to ya\nFirst we gonna rock, then we gonna roll\nThen we let it pop, go, let it go\nHES gon' give it to ya (uh), hes gon' give it to ya\nHES gon' give it to ya (uh), hes gon' give it to ya",
     "HES.....", "There HES goes. There HES goes again"]
    if 'hes' in checkMe.lower():

        print(f'\n\n{ran_ls[randint(0, len(ran_ls)-1)]}\n\n')


   

if __name__ == "__main__":

    import argparse
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("network_file")
    parser.add_argument("-r", type=str, default=False, help="To find all targets of a single regulator in the network, provide one regulator")
    parser.add_argument("-tt", type=str, default=False, help="Can be used with -r to filter targets by type \'miRNA\' or \'protein_coding\'")

    parser.add_argument("-t", type=str, default=False, help="To find all regulators of a single target in the network, provide one target")
    parser.add_argument("-rt", type=str, default=False, help="Can be used with -t to filter regulator by type \'miRNA\' or \'protein_coding\'")

    parser.add_argument("-g", type=str, default=False, help="find all targets and regulators of a single gene in the network, provide one gene")

    parser.add_argument("-cr", type=str, default=False, help="Requires ct (target), checks relationship between a target and regulator, will allow partial match")
    parser.add_argument("-ct", type=str, default=False, help="Requires cr (regulator), checks relationship between a target and regulator, will allow partial match")

    args = parser.parse_args()

    outpath = os.path.join(os.path.split(args.network_file)[0], 'catchOfTheDay')
    
    with open("fishNet.txt") as fishy:
        for line in fishy:
            print(line.strip())
    print("\nA Thomas Minchington script.\n")
    if not os.path.isdir(outpath):
        os.makedirs(outpath)
        exit()

    if args.r:
        
        if not args.tt:
            tt = 'all'
        else:
            tt= args.tt

        print('\n----------------------------------------------------\n')
        print(f"\n\nRunning in REGULATOR mode. Looking for {tt} targets of {args.r}. \n\n")
        random_snarky_comment(args.r)
        print('\n----------------------------------------------------\n')

        get_regulator_targets(args.network_file, [args.r], args.tt, outpath)

        print(f"All done! Data can be found in folder:\n{outpath}")
        exit()

    if args.t:

        if not args.rt:
            rt = 'all'
        else:
            rt= args.rt

        print('\n----------------------------------------------------\n')
        print(f"\n\nRunning in TARGET mode. Looking for {rt} regulators of {args.t}. \n\n")
        random_snarky_comment(args.t)
        print('\n----------------------------------------------------\n')

        get_target_regulators(args.network_file, [args.t], args.rt, outpath)

        print(f"All done! Data can be found in folder:\n{outpath}")
        exit()


    if args.g:

        print('\n----------------------------------------------------\n')
        print(f"\n\nRunning in GENE mode. Looking for {args.g} regulators and targets. \n\n")
        random_snarky_comment(args.g)
        print('\n----------------------------------------------------\n')

        surroundingNetwork(args.network_file, [args.g], outpath)

        print(f"All done! Data can be found in folder:\n{outpath}")
        exit()

    if args.cr and args.ct:

        print('\n----------------------------------------------------\n')
        print(f"\n\nRunning in pair check mode. Looking for interactions between {args.cr} and {args.ct}\n\n")
        random_snarky_comment(args.cr)
        random_snarky_comment(args.ct)
        print('\n----------------------------------------------------\n')

        checkCrCt(args.network_file, args.cr, args.ct, outpath)

        print(f"All done! Data can be found in folder:\n{outpath}")
        exit()

    elif args.cr or args.ct:
        exit("ERROR: cr and ct are both required to check a relationship")
    print('\n----------------------------------------------------\n')
    print("\n\nThat's a nice network you have there. Do you want to do something?\n\nTry running \'python fishingNet.py -h\' to see options\n\n")
    print('\n----------------------------------------------------\n')
    
