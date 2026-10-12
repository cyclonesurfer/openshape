# Bundled fonts

The six faces here are the Liberation fonts, version 2.1.5, unmodified:
Liberation Sans, Liberation Serif and Liberation Mono, each regular and
bold. They are licensed under the SIL Open Font License 1.1 (OFL.txt, with
the copyright notices); the OpenShape code that uses them is
GPL-3.0-or-later. ADR 0015 records why they were chosen.

Provenance: the files were taken from Debian's `fonts-liberation` package,
version 1:2.1.5-3, which builds them unchanged from the upstream sources at
https://github.com/liberationfonts/liberation-fonts (release 2.1.5).

| file | SHA-256 |
|---|---|
| LiberationSans-Regular.ttf | 4659bc0c58c5028dd488ec928d41d9265db43d9b669fc14ca8b0832daca7b144 |
| LiberationSans-Bold.ttf | 3973aa5054fb467dd5627245d3dc82e37bf16fe075756156a570455871351582 |
| LiberationSerif-Regular.ttf | 705903ae1382f0150f115d6f10a30f5b2f3d6ade649c474afe6366c39174fb65 |
| LiberationSerif-Bold.ttf | 9e66c25e20868756bf12eb24a520a4552727c9cd108f577f7ae33e3c0110e39c |
| LiberationMono-Regular.ttf | 395fa5ab8d40c8eba390ced528744ea75a7f69aabf3e68b6f925ca0e39a27370 |
| LiberationMono-Bold.ttf | 626655e94dd82f3f42549daf995c921b0915fa8ab1f4b839559e8892ea41d240 |

A test (`tests/fonts.rs`) checks these hashes, so the files cannot change
unnoticed. Documents store the outlines laid out from them (ADR 0015), so
replacing a font later would change only text laid out afterwards.
