# Inverse variance-stabilising transform

Inverts
[`.g_transform()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-g_transform.md):
`g_inverse(g_transform(I, c), c) == I`.

## Usage

``` r
.g_inverse(g, scale_c)
```

## Arguments

- g:

  Numeric vector on the transformed scale.

- scale_c:

  Single positive scale `c` (must match the value used in
  [`.g_transform()`](https://vorpalvorpal.github.io/hydroSense/reference/dot-g_transform.md)).

## Value

Numeric vector `scale_c * sinh(g)`.
