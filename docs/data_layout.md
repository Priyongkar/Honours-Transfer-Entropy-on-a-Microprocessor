# Data layout expected under HONOURS_DATA_ROOT

| Path | What it is | Used by |
|---|---|---|
| `data/raw/donkeykong.5000.ts.spikes.npy` | Published recording (639,999 x 3,510 uint8, change coding) | main run, 02, 04, 05 |
| `data/visual6502/` | Netlist (`transdefs.js`, `segdefs.js`, `nodenames.js`) | all |
| `data/regen/*.regen.ts.spikes.npy` | Recordings regenerated with Sim2600 (Donkey Kong, Pitfall, Space Invaders) | 06, 07 |
| `results/v5_bitcorrect_every_third_mi_ais_te_k1_12_L1_d1_ais_x1/` | Main run: groups, histories, score matrices, IDTxl networks, chip states and clamps | 02-05, 07, 08 |
| `vendor/Sim2600/` | Sim2600 (github.com/ericmjonas/Sim2600, commit f719877) with `simulator/` patch applied and built | 06, 08 |
