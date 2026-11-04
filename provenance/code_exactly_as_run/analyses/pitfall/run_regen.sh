#!/bin/bash
cd "/Users/priyongkarkundu/Documents/Honours 2"
caffeinate -i /Users/priyongkarkundu/opt/anaconda3/envs/honours6502/bin/python analyses/pitfall/record_game.py DonkeyKong 639999 data/regen/donkeykong.regen.ts.spikes.npy > analyses/pitfall/output/regen_dk.log 2>&1 &
caffeinate -i /Users/priyongkarkundu/opt/anaconda3/envs/honours6502/bin/python analyses/pitfall/record_game.py Pitfall 639999 data/regen/pitfall.regen.ts.spikes.npy > analyses/pitfall/output/regen_pitfall.log 2>&1 &
caffeinate -i /Users/priyongkarkundu/opt/anaconda3/envs/honours6502/bin/python analyses/pitfall/record_game.py SpaceInvaders 639999 data/regen/spaceinvaders.regen.ts.spikes.npy > analyses/pitfall/output/regen_si.log 2>&1 &
wait
