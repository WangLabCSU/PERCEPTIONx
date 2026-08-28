# Build a clone key from patient and clone id

The inverse of
[`parse_clone_keys()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/parse_clone_keys.md):
pastes patient and clone identifiers into a single `"patient@clone_id"`
key. The `"@"` separator is the canonical clone-key format used
throughout PERCEPTIONx.

## Usage

``` r
build_clone_key(patient, clone_id)
```

## Arguments

- patient:

  Character vector of patient IDs.

- clone_id:

  Character vector of clone IDs.

## Value

Character vector of `"patient@clone_id"` keys.

## See also

[`parse_clone_keys()`](https://wanglabcsu.github.io/PERCEPTIONx/reference/parse_clone_keys.md)
for splitting keys back apart.
