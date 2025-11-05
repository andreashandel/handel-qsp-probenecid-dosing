


## model-simulator-function.R
- **State re-assignment inside the virus clamp.** The block that forces `V` to zero for small values assigns `V <- 0` inside `with(...)` (lines 35–47). That local assignment does not feed back into the state vector unless the ODE solver revisits the callback immediately, so the following derivatives still see the pre-clamped `V`. If the intent is to clamp both the state and its derivative, consider returning the clamped value via the derivative (e.g., set `dV` based solely on the clamp) or manipulating `y` in the event function instead.
