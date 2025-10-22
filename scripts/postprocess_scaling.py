#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def is_nice_factorization(N, R):
    divs = [d for d in range(1, R+1) if R % d == 0]
    for px in divs:
        if N % px != 0: 
            continue
        for py in divs:
            if (R // px) % py != 0:
                continue
            pz = R // (px*py)
            if pz < 1:
                continue
            if N % py == 0 and N % pz == 0:
                return True
    return False

def main(path='build-scaling/results_scaling.csv', outdir='build-scaling/scaling_reports'):
    csv = Path(path)
    if not csv.exists():
        raise SystemExit(f'{csv} not found')
    out = Path(outdir); out.mkdir(exist_ok=True, parents=True)
    df = pd.read_csv(csv)

    for col in ['grid_n','ranks','omp_threads']:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    if 'wall_s' in df.columns:
        df['wall_s'] = pd.to_numeric(df['wall_s'], errors='coerce')
    if 'mean_us' in df.columns:
        df['mean_us'] = pd.to_numeric(df['mean_us'], errors='coerce')
    else:
        df['mean_us'] = np.nan
    if 'stddev_us' in df.columns:
        df['stddev_us'] = pd.to_numeric(df['stddev_us'], errors='coerce')
    else:
        df['stddev_us'] = np.nan

    # prefer mean_us, else fallback to wall_s
    if df['mean_us'].isna().all() and 'wall_s' in df.columns:
        df['mean_us'] = df['wall_s'] * 1e6

    df['nice_factorization'] = [is_nice_factorization(int(n), int(r)) if pd.notna(n) and pd.notna(r) else False
                                for n, r in zip(df.get('grid_n', []), df.get('ranks', []))]

    agg = (df.groupby(['grid_n','ranks'], dropna=True)
             .agg(mean_us=('mean_us','mean'),
                  n=('mean_us','count'),
                  nice=('nice_factorization','max'))
             .reset_index()
             .sort_values(['grid_n','ranks']))

    rows = []
    for g, gdf in agg.groupby('grid_n'):
        gdf = gdf.sort_values('ranks')
        base_time = gdf['mean_us'].iloc[0]
        base_ranks = gdf['ranks'].iloc[0]
        for _, row in gdf.iterrows():
            r = row['ranks']; t = row['mean_us']
            speed = base_time / t if t and t>0 else np.nan
            eff = speed / (r / base_ranks) if r and base_ranks else np.nan
            rows.append(dict(grid_n=g, ranks=r, mean_us=t,
                             speedup_vs_min_ranks=speed,
                             efficiency_vs_min_ranks=eff,
                             nice_factorization=bool(row['nice']),
                             samples=int(row['n'])))
    res = pd.DataFrame(rows)
    cleaned = out / 'results_scaling_cleaned.csv'
    res.to_csv(cleaned, index=False)

    for g, gdf in res.groupby('grid_n'):
        gdf = gdf.sort_values('ranks')
        plt.figure()
        plt.plot(gdf['ranks'], gdf['speedup_vs_min_ranks'], marker='o')
        plt.xscale('log', base=2)
        plt.xlabel('MPI ranks'); plt.ylabel('Speedup (vs min ranks)')
        plt.title(f'Strong scaling speedup — N={int(g)}^3')
        plt.grid(True, which='both', linestyle='--', alpha=0.5)
        plt.savefig(out / f'speedup_N{int(g)}.png', bbox_inches='tight')
        plt.close()

        plt.figure()
        plt.plot(gdf['ranks'], gdf['efficiency_vs_min_ranks'], marker='o')
        plt.xscale('log', base=2)
        plt.xlabel('MPI ranks'); plt.ylabel('Parallel efficiency')
        plt.title(f'Parallel efficiency — N={int(g)}^3')
        plt.grid(True, which='both', linestyle='--', alpha=0.5)
        plt.savefig(out / f'efficiency_N{int(g)}.png', bbox_inches='tight')
        plt.close()

    print(f'Wrote cleaned CSV: {cleaned}')
    print(f'Plots in: {out}')

if __name__ == '__main__':
    main()
