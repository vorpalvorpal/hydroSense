# .g_transform() / .g_inverse(): errors on a non-positive scale c

    Code
      leachatetools:::.g_transform(1, 0)
    Condition
      Error in `.assert_pos_scale()`:
      ! `scale_c` must be a single positive number.

---

    Code
      leachatetools:::.g_transform(1, -3)
    Condition
      Error in `.assert_pos_scale()`:
      ! `scale_c` must be a single positive number.

