# CLAUDE.md

Instructions for Claude Code when working with this repository.

## Project Overview

CHOWTape is an audio effect module for Move Anything that provides analog tape modeling based on the Jiles-Atherton hysteresis model. It is a pure C port of the [AnalogTapeModel](https://github.com/jatinchowdhury18/AnalogTapeModel) by Jatin Chowdhury.

**License: GPL-3.0** (derived from GPL-3.0 original)

## Architecture

```
src/
  dsp/
    chowtape.c          # Main DSP implementation (all stages)
    audio_fx_api_v1.h   # Audio FX API (from move-anything)
    plugin_api_v1.h     # Plugin API types (from move-anything)
  module.json           # Module metadata
```

## DSP Chain

1. **Tone Control** - Pre-EQ with low/high shelving filters at 500Hz crossover
2. **Hysteresis** - Jiles-Atherton magnetic hysteresis model with RK4 solver (the core tape saturation)
3. **DC Blocker** - Removes DC offset introduced by hysteresis
4. **Loss Filter** - Head loss modeling (1-pole lowpass + head bump resonance), varies with tape speed
5. **Chew** - Tape damage: random dropout events (power distortion) + lowpass filtering
6. **Degrade** - Tape degradation: noise injection + lowpass + gain reduction
7. **Wow & Flutter** - LFO-modulated delay line with Ornstein-Uhlenbeck process for variance
8. **Dry/Wet Mix** - Latency-compensated parallel dry path

## Parameters

| Parameter | Range | Default | Description |
|-----------|-------|---------|-------------|
| drive | 0-1 | 0.5 | Hysteresis input drive |
| saturation | 0-1 | 0.5 | Magnetic saturation |
| bias | 0-1 | 0.5 | Tape bias / width |
| tone | 0-1 | 0.5 | Tilt EQ (dark to bright) |
| speed | 0-1 | 0.5 | Tape speed (affects HF loss) |
| flutter | 0-1 | 0.0 | Combined wow & flutter depth |
| chew | 0-1 | 0.0 | Tape chew/damage amount |
| degrade | 0-1 | 0.0 | Degradation (noise + filter) |
| mix | 0-1 | 1.0 | Dry/wet balance |
| output | 0-1 | 0.5 | Output level |

## Key Implementation Details

### Hysteresis Model
- Jiles-Atherton model parameters: M_s (saturation), a (shape), k (pinning), c (reversibility), alpha (coupling)
- cook() maps drive/width/saturation to physical parameters
- RK4 (Runge-Kutta 4th order) numerical solver
- Alpha derivative approximation (weighted 0.75 coefficient)
- Double precision for numerical stability

### Build Commands

```bash
./scripts/build.sh      # Build for ARM64 via Docker
./scripts/install.sh    # Deploy to Move
```

## Signal Chain Integration

Module declares `"chainable": true` and `"component_type": "audio_fx"` in module.json.
Installs to: `/data/UserData/move-anything/modules/audio_fx/chowtape/`
