# CHOWTape for Schwung

Analog tape model audio effect for [Schwung](https://github.com/charlesvestal/schwung), based on the [AnalogTapeModel](https://github.com/jatinchowdhury18/AnalogTapeModel) by Jatin Chowdhury.

The core of the effect is a Jiles-Atherton magnetic hysteresis model, solved with a 4th-order Runge-Kutta integrator. This physically models how magnetic tape saturates, producing warm harmonic distortion and natural compression that varies with input level.

## Features

- **Hysteresis saturation** — Jiles-Atherton model with RK4 solver for physically accurate tape distortion
- **Tape speed / loss filter** — Head loss modeling with speed-dependent HF rolloff and head bump resonance
- **Wow & flutter** — LFO-modulated delay with Ornstein-Uhlenbeck process for realistic pitch drift
- **Chew** — Tape damage simulation with random dropout events and filtering
- **Degradation** — Noise injection and lowpass filtering for worn tape character
- **Tone control** — Pre-tape tilt EQ with low and high shelving filters

## Parameters

| Knob | Description |
|------|-------------|
| Drive | Input level into the hysteresis model |
| Saturation | Magnetic saturation — lower = thicker/more compressed |
| Bias | Tape bias / width — changes harmonic character |
| Tone | Pre-tape tilt EQ (dark to bright) |
| Speed | Tape speed — affects HF loss and head bump |
| Flutter | Combined wow & flutter depth |
| Chew | Tape damage / dropout amount |
| Mix | Dry/wet balance |

Additional parameters (accessible via chain editor): Degrade, Output.

## Build

```bash
./scripts/build.sh      # Cross-compile for ARM64 via Docker
./scripts/install.sh    # Deploy to Move via SSH
```

## Credits

- **AnalogTapeModel** by [Jatin Chowdhury](https://github.com/jatinchowdhury18) — original DSP algorithms
- **Hysteresis model** from [DAFx paper](https://ccrma.stanford.edu/~jatin/420/tape/TapeModel_DAFx.pdf) by Jatin Chowdhury (Stanford CCRMA)

## License

GPL-3.0-or-later — see [LICENSE](LICENSE) for details.

This is a derivative work of AnalogTapeModel, which is licensed under GPL-3.0.
