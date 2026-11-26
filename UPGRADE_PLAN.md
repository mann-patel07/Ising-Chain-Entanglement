# Plan for periodic chains with free fermions

Date: 2026-11-26

Written before these results existed. The pre-registered estimate (open chains, L = 1024 and 2048) and its test are unchanged.

## Question
Do large periodic chains give c = 1/2, and how fast do they converge compared with open chains?

## Method
- In the even parity sector, which holds the ground state for h > 0, the periodic spin chain maps to Majorana fermions with an antiperiodic boundary: the bond between the last and the first site enters with the opposite sign.
- Check against exact diagonalization of periodic chains (L = 4 to 12, h = 0.5, 1, 2): energy to 1e-9 and every block entropy to 1e-8.
- At h = 1: half-chain entropy for L = 16, 32, ..., 2048, and c from each doubling, c = 3 (S(2L) - S(L)) / ln 2.
- A fit of S(l) against (c/3) ln[(L/pi) sin(pi l / L)] at L = 1024, l = 64 to 512 in steps of 16.
- Supporting estimates. The pre-registered main estimate stays the main estimate. I report the results whatever they are.
