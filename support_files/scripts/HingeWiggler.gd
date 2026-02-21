@tool
class_name HingeWiggler extends Node3D

# without the factor we would have to set compliance in the range of
# 0.00001 to 0.001...
const COMPLIANCE_FACTOR : float = 0.00001

@export var mass: float = 0.5:
	set(value): 
		if _tip_particle != -1:
			WiggleKit.particle_set_mass(_tip_particle, value)
		mass = value

@export var damping : float = 0.5:
	set(value): 
		if _tip_particle != -1:
			WiggleKit.particle_set_damping(_tip_particle, value)
		damping = value

@export var width: float = 0.5
@export var height: float = 1.0

@export var compliance : float = 0:
	set(value):
		for dist_constraint in _distance_constraints:
			WiggleKit.distance_constraint_set_compliance(dist_constraint, value)
		compliance = value

@export var wiggle_intensity : float = 1

@export var gravity : Vector3 = Vector3(0,-10,0):
	set(value): 
		if _tip_particle != -1:
			WiggleKit.particle_set_gravity(_tip_particle, value)
		gravity = value

@export var Colliders: Array[WiggleShapeBase] = []
@export var ShapeDistanceConstraints: Array[WiggleShapeBase] = []

var _position_particle1 : int = -1
var _position_particle2 : int = -1
var _tip_particle : int = -1

var _distance_constraints: Array[int] = []

func _enter_tree() -> void:
	if Engine.is_editor_hint():
		return
	
	# position particles will be kept statically at our global_position (base of hinge)
	_position_particle1 = WiggleKit.particle_create(to_global(Vector3.LEFT * width), 0, Vector3.ZERO, 0)
	_position_particle2 = WiggleKit.particle_create(to_global(Vector3.RIGHT * width), 0, Vector3.ZERO, 0)
	
	# tip particle
	_tip_particle = WiggleKit.particle_create(to_global(Vector3.DOWN * height), mass, gravity, damping)
	
	# distance constraints from position particles to _tip particle
	_distance_constraints.append(WiggleKit.distance_constraint_create(_position_particle1, _tip_particle, 0))
	_distance_constraints.append(WiggleKit.distance_constraint_create(_position_particle2, _tip_particle, 0))
	
	for collider in Colliders:
		if collider: collider.add_particle_to_collider(_tip_particle)
	for sdf in ShapeDistanceConstraints:
		if sdf: sdf.add_particle_to_sdf_constraint(_tip_particle)

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if WiggleKit:
		for c_id in _distance_constraints:
			WiggleKit.distance_constraint_free(c_id)
		_distance_constraints.clear()
		
		if _position_particle1 != -1: WiggleKit.particle_free(_position_particle1)
		if _position_particle2 != -1: WiggleKit.particle_free(_position_particle2)
		if _tip_particle != -1: WiggleKit.particle_free(_tip_particle)
		
		_position_particle1 = -1
		_position_particle2 = -1
		_tip_particle = -1

func _physics_process(delta: float) -> void:
	if Engine.is_editor_hint():
		return
	
	if _tip_particle == -1:
		return
	
	var tip_pos := WiggleKit.particle_get_position(_tip_particle)
	
	# Y axis points from global_position to the tip
	var y_axis := (global_position - tip_pos).normalized()
	# a bit of lerping, so that the global changes are not completely overwritten
	y_axis = global_basis.y.slerp(y_axis, 1 - exp(-wiggle_intensity * delta))
	var x_axis := global_basis.x
	var z_axis := x_axis.cross(y_axis).normalized()
	global_basis = Basis(x_axis, y_axis, z_axis).orthonormalized()
	
	# Update the static position particles
	WiggleKit.particle_set_position(_position_particle1, to_global(Vector3.LEFT * width))
	WiggleKit.particle_set_position(_position_particle2, to_global(Vector3.RIGHT * width))
	
	# project the tip_pos into the rotation plane (otherwise there might be huge velocities)
	var local_tip_pos := to_local(tip_pos)
	local_tip_pos.x = 0
	tip_pos = to_global(local_tip_pos)
	WiggleKit.particle_set_position(_tip_particle, tip_pos)
