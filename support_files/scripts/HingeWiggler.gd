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

@export_range(0, 1, 0.05) var ParentMovementInfluence : float = 0.9
## Low-pass factor for the parent-motion estimate (1 = no smoothing).
@export_range(0.01, 1, 0.01) var MovementSmoothing : float = 0.4
## Clamps the per-tick parent displacement in metres to reject network snaps (0 = off).
@export var MaxMovementPerTick : float = 0.5

@export var width: float = 0.5
@export var height: float = 1.0

@export var compliance : float = 0:
	set(value):
		for dist_constraint in _distance_constraints:
			WiggleKit.distance_constraint_set_compliance(dist_constraint, value)
		compliance = value


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

var _parent_node3D : Node3D
var _parent_motion := WiggleParentMotion.new()

func _enter_tree() -> void:
	if Engine.is_editor_hint():
		return

	add_to_group("wigglers")
	_parent_node3D = get_parent_node_3d()

	# position particles will be kept statically at our position (base of hinge, parent-local)
	_position_particle1 = WiggleKit.particle_create(transform * (Vector3.LEFT * width), 999, Vector3.ZERO, 0)
	_position_particle2 = WiggleKit.particle_create(transform * (Vector3.RIGHT * width), 999, Vector3.ZERO, 0)

	# tip particle
	_tip_particle = WiggleKit.particle_create(transform * (Vector3.DOWN * height), mass, gravity, damping)

	# distance constraints from position particles to _tip particle
	_distance_constraints.append(WiggleKit.distance_constraint_create(_position_particle1, _tip_particle, 0))
	_distance_constraints.append(WiggleKit.distance_constraint_create(_position_particle2, _tip_particle, 0))

	for collider in Colliders:
		if collider: collider.add_particle_to_collider(_tip_particle)
	for sdf in ShapeDistanceConstraints:
		if sdf: sdf.add_particle_to_sdf_constraint(_tip_particle)

	_parent_motion.reset(_parent_node3D)

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

## Drops parent-motion history; call after a hard parent teleport/snap.
func reset_motion() -> void:
	if _parent_node3D:
		_parent_motion.reset(_parent_node3D)

func _physics_process(_delta: float) -> void:
	if Engine.is_editor_hint():
		return

	if _tip_particle == -1:
		return

	# Update the static position particles
	WiggleKit.particle_set_position(_position_particle1, transform * (Vector3.LEFT * width))
	WiggleKit.particle_set_position(_position_particle2, transform * (Vector3.RIGHT * width))

	var to_local_rotation := _parent_node3D.global_basis.inverse()
	var parent_accel := _parent_motion.update(_parent_node3D, to_local_rotation, MovementSmoothing, MaxMovementPerTick)

	WiggleKit.particle_set_gravity(_tip_particle, to_local_rotation * gravity)
	var tip_pos := WiggleKit.particle_get_position(_tip_particle)
	if ParentMovementInfluence > 0:
		tip_pos -= parent_accel * ParentMovementInfluence

	# Y axis points from our (parent-local) position to the tip
	var y_axis := (position - tip_pos).normalized()
	var x_axis := basis.x
	var z_axis := x_axis.cross(y_axis).normalized()
	basis = Basis(x_axis, y_axis, z_axis).orthonormalized()

	# project the tip_pos into the rotation plane (otherwise there might be huge velocities)
	var to_local_transform := transform.affine_inverse()
	var local_tip_pos := to_local_transform * tip_pos
	local_tip_pos.x = 0
	tip_pos = transform * local_tip_pos
	WiggleKit.particle_set_position(_tip_particle, tip_pos)
