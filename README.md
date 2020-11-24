-------------
fishingNET.py
--------------

Used to look at specific interactions within the network from Minchington et al. 2020. Dynamical gene regulatory networks are tuned by transcriptional autoregulation with microRNA feedback. Sci Rep. 2020 Jul 31;10(1):12960. doi: 10.1038/s41598-020-69791-5. PMID: 32737375; PMCID: PMC7395740.

usage: fishingNet.py [-h] [-r R] [-tt TT] [-t T] [-rt RT] [-g G] [-cr CR]
                     [-ct CT]
                     network_file

positional arguments:
  network_file

optional arguments:
  -h, --help    show this help message and exit
  -r R          To find all targets of a single regulator in the network,
                provide one regulator
  -tt TT        Can be used with -r to filter targets by type 'miRNA' or
                'protein_coding'
  -t T          To find all regulators of a single target in the network,
                provide one target
  -rt RT        Can be used with -t to filter regulator by type 'miRNA' or
                'protein_coding'
  -g G          find all targets and regulators of a single gene in the
                network, provide one gene
  -cr CR        Requires ct (target), checks relationship between a target and
                regulator, will allow partial match
  -ct CT        Requires cr (regulator), checks relationship between a target
                and regulator, will allow partial match
