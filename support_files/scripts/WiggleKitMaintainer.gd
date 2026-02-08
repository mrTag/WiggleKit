extends Node

var _last_simulation_duration_ns : float

func _ready() -> void:
	# some of the wiggle depends on the position of other physics objects,
	# so we take care to do our calculations last
	process_physics_priority = 999
	Performance.add_custom_monitor("WiggleKit step duration ns", func(): return _last_simulation_duration_ns)

func _physics_process(delta: float) -> void:
	var start_time := Time.get_ticks_usec()
	WiggleKit.iteration(delta)
	_last_simulation_duration_ns = Time.get_ticks_usec() - start_time
