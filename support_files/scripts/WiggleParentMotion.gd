class_name WiggleParentMotion extends RefCounted

## Tracks a wiggler parent's motion and produces the local-space *acceleration* term used
## to inject inertia into the simulation.
##
## The wiggler particles live in the parent's local space, so a parent moving at a constant
## velocity should NOT deflect them (a pendulum at constant speed hangs straight down) -
## only a *change* of velocity (acceleration) should. We therefore return the change in the
## per-tick local displacement between ticks (≈ a·dt²), not the displacement itself (≈ v·dt).
##
## The estimate is low-pass filtered and clamped so a sudden parent position jump (e.g. a
## multiplayer network correction) doesn't inject a violent spike. On a genuine teleport,
## call reset() to drop the history entirely.

var _previous_global_pos : Vector3
var _smoothed_movement : Vector3      # low-pass filtered local displacement (≈ v·dt)
var _initialized : bool = false


## Drop all motion history and re-anchor to the parent's current position. Call on first
## use and whenever the parent teleports / snaps so no fake acceleration is injected.
func reset(parent: Node3D) -> void:
	_previous_global_pos = parent.global_position
	_smoothed_movement = Vector3.ZERO
	_initialized = true


## Returns the filtered local-space acceleration term to subtract from particle positions.
## `to_local_rotation` is the parent's inverse global basis (the caller already computes it
## for gravity, so it's passed in to avoid a second inverse).
## `smoothing` in (0,1] is the low-pass factor (1 = no smoothing); `max_step` in metres
## clamps the per-tick displacement (<= 0 disables the clamp).
func update(parent: Node3D, to_local_rotation: Basis, smoothing: float, max_step: float) -> Vector3:
	if not _initialized:
		reset(parent)
		return Vector3.ZERO

	var raw := to_local_rotation * (parent.global_position - _previous_global_pos)  # ≈ v·dt (local)
	_previous_global_pos = parent.global_position

	if max_step > 0.0 and raw.length() > max_step:        # reject network snaps
		raw = raw.normalized() * max_step

	var prev := _smoothed_movement
	_smoothed_movement = _smoothed_movement.lerp(raw, smoothing)   # low-pass the velocity estimate
	return _smoothed_movement - prev                       # ≈ a·dt² (filtered)
