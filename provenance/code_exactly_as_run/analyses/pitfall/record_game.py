"""Regenerate a Jonas & Kording (2017) transistor recording with their simulator (Sim2600, CC0),
ported to Python 3 in vendor/Sim2600 (lib2to3; pickle latin1; relative imports; ROM bytes; build with
setup_py3.py). Recording format, matched exactly to data/raw/donkeykong.5000.ts.spikes.npy:
row i = 1 where a transistor's state after TIA half-clock i+1 differs from its state after half-clock i
(change coding), shape [rows, 3510] uint8, columns in the simulator's transistor order.

Check (4 Nov 2026): the first 2,990 rows regenerated for Donkey Kong match the published file in every
one of 10,494,900 entries; turn-on coding gives 386,469 mismatches.

Usage: python record_game.py GAME ROWS OUT.npy      GAME in {DonkeyKong, Pitfall, SpaceInvaders}
"""
import sys, time, hashlib, json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'vendor/Sim2600'))
from sim2600 import params, sim6502, simTIA            # noqa: E402
from sim2600.sim2600Console import Sim2600Console       # noqa: E402

ROMS = {'DonkeyKong': params.ROMS_DONKEY_KONG, 'Pitfall': params.ROMS_PITFALL,
        'SpaceInvaders': params.ROMS_SPACE_INVADERS}


def main(game, rows, out):
    rom = ROMS[game]
    c = Sim2600Console(rom, sim6502.MySim6502, simTIA.MySimTIA)
    cpu = c.sim6507
    n = len(cpu.getTransistorState())
    mm = np.lib.format.open_memmap(out, mode='w+', dtype=np.uint8, shape=(rows, n))
    c.advanceOneHalfClock(); prev = cpu.getTransistorState().copy()
    t0 = time.time()
    for i in range(rows):
        c.advanceOneHalfClock(); cur = cpu.getTransistorState()
        mm[i] = (cur != prev); prev = cur.copy()
        if (i + 1) % 50000 == 0:
            print(f'{game}: {i + 1}/{rows} rows, {time.time() - t0:.0f} s', flush=True)
    mm.flush(); del mm
    meta = dict(game=game, rom=Path(rom).name, rom_sha256=hashlib.sha256(Path(rom).read_bytes()).hexdigest(),
                rows=rows, transistors=n, seconds=round(time.time() - t0, 1))
    json.dump(meta, open(str(out) + '.json', 'w'), indent=1)
    print('done', meta, flush=True)


if __name__ == '__main__':
    import builtins, io, contextlib
    game, rows, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    _print = builtins.print
    builtins.print = lambda *a, **k: _print(*a, **k) if str(a[0]).startswith((game, 'done')) else None
    main(game, rows, out)
