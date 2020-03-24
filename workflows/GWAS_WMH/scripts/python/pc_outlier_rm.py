# code from Ken Hanscombe

import numpy as np
import pandas as pd
import subprocess
import argparse
import sys

parser = argparse.ArgumentParser(
    description='Iterative outlier removal with PLINK 2.0.')

parser.add_argument('--plink2', default=None,
                    type=str, help='Path to PLINK 2.0.')
parser.add_argument('--bfile', default=None, type=str,
                    help='PLINK binary file set prefix.')
parser.add_argument('--out', default='plink', type=str,
                    help='Output file path (default: plink).')
parser.add_argument('--npc-rm', default=10, type=int,
                    help='Number of PCs to remove outliers on (default: 10).')
parser.add_argument('--nsd-rm', default=6, type=int,
                    help='Number of SDs from mean define outliers (default: 6).')


def _outliers(pc_series, nsd):
    """Identifies genetic principal component outliers.
    """
    pc = pc_series
    sd = np.std(pc)
    mu = np.mean(pc)
    out = (mu - pc).abs() > (nsd * sd)
    nout = out.sum()
    print(pc.name, ':', nout, 'outliers')
    return list(pc_series.index.values[out])


def _outlier_ids(bfile, plink2, nsd, npc, out):
    """Runs PLINK 2.0 PCA and lists outliers along first 10 PCs.
    """
    print(f'\nRemoving outliers along first {npc} PCs.')
    print(f'Outliers defined as {nsd} SDs from the mean.')
    print('Calculating PCs...\n')

    try:
        subprocess.run(f'{plink2} --bfile {bfile} --pca approx --out {out}',
                       shell=True, check=True)
    except:
        sys.exit('\n\n!! PLINK did not exit cleanly. Possible system resource allocation error. Check output\n\n')

    pcs = pd.read_csv(f'{out}.eigenvec', sep=r'\s+', header=0)
    out_ids = [_outliers(pcs[p], nsd=nsd) for p in pcs.drop(
        ['#FID', 'IID'], axis=1).columns]
    out_ids = set(np.sort(np.concatenate(out_ids).ravel()))
    return out_ids


if __name__ == '__main__':
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()

    outlier_ids = _outlier_ids(f'{args.bfile}', plink2=f'{args.plink2}',
                               nsd=args.nsd_rm, npc=args.npc_rm,
                               out=f'{args.out}')

    if len(outlier_ids) == 0:
        print('\nNo PC outliers to remove.\n')
        bfile = args.bfile
        subprocess.run(f'{args.plink2} --bfile {bfile} --make-bed \
                       --out {args.out}', shell=True)

        print('\nDONE.')
        print(f'{args.out}.bed, {args.out}.bim, {args.out}.fam is the sample without outliers.')

    else:
        iteration = 1
        while len(outlier_ids) > 0:
            print(len(outlier_ids), 'outliers to remove.')
            print('Removing outliers...\n')
            pcs = pd.read_csv(f'{args.out}.eigenvec', sep=r'\s+', header=0)
            (pcs[['#FID', 'IID']].iloc[list(outlier_ids)]
                .to_csv(f'{args.out}_iter{iteration}.pc.out',
                        sep=' ', header=None, index=False))

            if iteration == 1:
                bfile = args.bfile
            else:
                bfile = args.out

            subprocess.run(f'{args.plink2} --bfile {bfile} \
                           --remove {args.out}_iter{iteration}.pc.out \
                           --make-bed --out {args.out}', shell=True)

            outlier_ids = _outlier_ids(f'{args.out}', plink2=f'{args.plink2}',
                                       nsd=args.nsd_rm, npc=args.npc_rm,
                                       out=f'{args.out}')
            iteration += 1

        print('\nDONE.')
        print(f'Fileset without outliers: {args.out}.bed, {args.out}.bim, {args.out}.fam')
        print(f'Samples removed at each iteration: {args.out}_iter<iteration>.pc.out\n')
