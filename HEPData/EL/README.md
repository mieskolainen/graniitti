# Elastic HEPData

Elastic-scattering HEPData JSON tables.

| Energy | Folder | Record | Process | Relevant content |
|---:|---|---:|---|---|
| 23.5 GeV | `HEPData-ins214689-v1-json` | `ins214689` | ISR `pp` | `Table1.json` |
| 53 GeV | `HEPData-ins212895-v1-json` | [ins212895](https://www.hepdata.net/record/ins212895) | ISR `pp`, `pbar p` | `Table1.json` (`pp`), `Table2.json` (`pbar p`), diffraction dip |
| 62.5 GeV | `HEPData-ins84176-v1-json` | `ins84176` | ISR `pp` | ISR low-energy `d sigma/dt` tables |
| 200 GeV | `HEPData-ins1791591-v1-json` | `ins1791591` | STAR `pp` | `Figure5.json` |
| 546 GeV | `HEPData-ins201990-v1-json` | `ins201990` | UA4/SPS `pbar p` | `Table3.json` |
| 546 GeV | `HEPData-ins359411-v1-json` | `ins359411` | CDF/E741 `pbar p` | `Table3.json` |
| 1.8 TeV | `HEPData-ins359411-v1-json` | `ins359411` | CDF/E741 `pbar p` | `Table4.json` |
| 1.96 TeV | `HEPData-ins1117021-v1-json` | `ins1117021` | D0 `pbar p` | `Table1.json` |
| 7 TeV | `HEPData-ins1220862-v1-json` | `ins1220862` | TOTEM `pp` | `Table1.json`, low `|t|` |
| 7 TeV | `HEPData-ins922651-v1-json` | `ins922651` | TOTEM `pp` | `Table1.json`, high `|t|` |
| 8 TeV | `HEPData-ins1489188-v1-json` | `ins1489188` | TOTEM `pp` | `Table1.json`, CNI |
| 13 TeV | `HEPData-ins1710340-v1-json` | `ins1710340` | TOTEM `pp` | `Table1.json` |

The original `ins212895` JSON tables are unchanged. The [ISR dip icepack](../../icepack/ELASTIC/ISR_212895/dataset.json) uses the published bin edges and overlays explicit `single` and `double` samples for each beam system. Its reader follows [Table I and p. 2181](https://doi.org/10.1103/PhysRevLett.54.2180): the point errors include statistics and correction uncertainties, the absolute normalization errors are 20% for `pp` and 30% for `pbar p`, and their relative normalization error is 20%. The relative error fixes the correlation between the two normalization scales. Run `bash icepack/run.sh ELASTIC/ISR_212895`.
