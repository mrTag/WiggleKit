@tool
class_name RopeWiggler extends MultiMeshInstance3D

const COMPLIANCE_FACTOR : float = 0.00001

@export var segment_count: int = 10:
	set(value):
		segment_count = max(1, value)

@export var total_length: float = 2.0

@export var segment_mass: float = 0.1:
	set(value):
		segment_mass = value
		var last_idx := segment_count - 1
		for i in _particles.size():
			var p_id : int = _particles[i]
			if p_id == -1: continue
			if end_attachment and i == last_idx: continue
			WiggleKit.particle_set_mass(p_id, value)

@export var segment_damping: float = 0.5:
	set(value):
		segment_damping = value
		for p_id in _particles:
			if p_id != -1: WiggleKit.particle_set_damping(p_id, value)

@export_range(0, 1, 0.05) var ParentMovementInfluence : float = 0.9
## Low-pass factor for the parent-motion estimate (1 = no smoothing).
@export_range(0.01, 1, 0.01) var MovementSmoothing : float = 0.4
## Clamps the per-tick parent displacement in metres to reject network snaps (0 = off).
@export var MaxMovementPerTick : float = 0.5

@export var gravity: Vector3 = Vector3(0, -10, 0):
	set(value):
		gravity = value
		for p_id in _particles:
			if p_id != -1: WiggleKit.particle_set_gravity(p_id, value)

@export var stretch_compliance: float = 0.0:
	set(value):
		stretch_compliance = value
		for dist_constraint in _distance_constraints:
			WiggleKit.distance_constraint_set_compliance(dist_constraint, value * COMPLIANCE_FACTOR)

@export var Colliders: Array[WiggleShapeBase] = []
@export var ShapeDistanceConstraints: Array[WiggleShapeBase] = []

@export var end_attachment: Node3D = null

# Fixed anchor particle (mass 0)
var _anchor_particle: int = -1

# Rope segment particles
var _particles: Array[int] = []

var _distance_constraints: Array[int] = []

var _parent_motion := WiggleParentMotion.new()
var _parent_node3D : Node3D

var init : bool = false
func _enter_tree() -> void:
	if Engine.is_editor_hint():
		_setup_multimesh()
		return

	add_to_group("wigglers")
	_parent_node3D = get_parent_node_3d()

	_setup_multimesh()
	_create_particles()
	_create_constraints()
	_register_colliders()
	_parent_motion.reset(_parent_node3D)
	
	# wait with the update until we are fully initialized
	init = false
	set_physics_process(false)
	await get_tree().process_frame
	if end_attachment:
		WiggleKit.particle_set_position(_particles[segment_count - 1], _parent_node3D.to_local(end_attachment.global_position))
		WiggleKit.warmup(200, _particles)
	
	set_physics_process(true)
	init = true


func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if WiggleKit:
		for c_id in _distance_constraints:
			if c_id != -1: WiggleKit.distance_constraint_free(c_id)
		_distance_constraints.clear() 
		
		if _anchor_particle != -1: WiggleKit.particle_free(_anchor_particle)
		
		for p_id in _particles:
			if p_id != -1: WiggleKit.particle_free(p_id)
		_particles.clear()


func _setup_multimesh() -> void:
	if multimesh == null:
		multimesh = MultiMesh.new()
		multimesh.transform_format = MultiMesh.TRANSFORM_3D
	multimesh.instance_count = segment_count


func _create_particles() -> void:
	var segment_length := total_length / float(segment_count)
	
	# Create anchor (fixed) particle at origin
	_anchor_particle = WiggleKit.particle_create(position, 0, Vector3.ZERO, 0)
	
	# Create rope segment particles
	for i in range(segment_count):
		var y_offset := -segment_length * (i + 1)
		var pos := transform * Vector3(0, y_offset, 0)
		var mass := segment_mass
		if end_attachment and i == segment_count - 1: mass = 0.0
		var p_id = WiggleKit.particle_create(pos, mass, gravity, segment_damping)
		_particles.append(p_id)


func _create_constraints() -> void:
	var stretch_comp := stretch_compliance * COMPLIANCE_FACTOR
	
	# Connect first particle to anchor
	_distance_constraints.append(WiggleKit.distance_constraint_create(
		_anchor_particle, _particles[0], stretch_comp))
	
	# Connect each particle to the next
	for i in range(segment_count - 1):
		_distance_constraints.append(WiggleKit.distance_constraint_create(
			_particles[i], _particles[i + 1], stretch_comp))


func _register_colliders() -> void:
	for p_id in _particles:
		for collider in Colliders:
			if collider: collider.add_particle_to_collider(p_id)
		for sdf in ShapeDistanceConstraints:
			if sdf: sdf.add_particle_to_sdf_constraint(p_id)


## Drops parent-motion history; call after a hard parent teleport/snap.
func reset_motion() -> void:
	if _parent_node3D:
		_parent_motion.reset(_parent_node3D)

func _physics_process(_delta: float) -> void:
	if not init:
		return
	if Engine.is_editor_hint():
		_update_editor_preview()
		return
	
	if _particles.size() < segment_count:
		return
	
	# Update anchor position to follow node
	WiggleKit.particle_set_position(_anchor_particle, position)
	# Update end attachment if set
	if end_attachment and is_instance_valid(end_attachment):
		var last_idx := segment_count - 1
		WiggleKit.particle_set_position(_particles[last_idx], _parent_node3D.to_local(end_attachment.global_position))
	
	var to_local_rotation := _parent_node3D.global_basis.inverse()
	var parent_accel := _parent_motion.update(_parent_node3D, to_local_rotation, MovementSmoothing, MaxMovementPerTick)
	var dynamic_particle_count := _particles.size() if not end_attachment else _particles.size()-1
	for i in dynamic_particle_count:
		var p_id : int = _particles[i]
		WiggleKit.particle_set_gravity(p_id, to_local_rotation * gravity)
		if ParentMovementInfluence > 0:
			var particle_pos := WiggleKit.particle_get_position(p_id)
			WiggleKit.particle_set_position(p_id, particle_pos - parent_accel * ParentMovementInfluence)
	

func _process(_delta:float) -> void:
	# Update multimesh instances
	_update_multimesh()


func _update_multimesh() -> void:
	if multimesh == null or _particles.size() < segment_count:
		return
	
	var to_local_transform := transform.affine_inverse()
	var prev_pos: Vector3 = to_local_transform * WiggleKit.particle_get_position(_anchor_particle) 
	for i in range(segment_count):
		var curr_pos := to_local_transform * WiggleKit.particle_get_position(_particles[i])
		
		# Calculate basis from direction
		var y_axis := (prev_pos - curr_pos).normalized()
		
		if y_axis.length_squared() < 0.0001:
			y_axis = Vector3.UP
		
		# Create orthonormal basis
		var x_axis := Vector3.RIGHT if abs(y_axis.x) < 0.9 else Vector3.FORWARD
		x_axis = (x_axis - y_axis * x_axis.dot(y_axis)).normalized()
		var z_axis := x_axis.cross(y_axis).normalized()
		
		var basis := Basis(x_axis, y_axis, z_axis).orthonormalized()
		
		# Position at midpoint between this and previous particle
		var mid_pos := (prev_pos + curr_pos) / 2.0
		
		var xform := Transform3D(basis, mid_pos)
		multimesh.set_instance_transform(i, xform)
		
		prev_pos = curr_pos
		


func _update_editor_preview() -> void:
	if multimesh == null:
		return
	
	var segment_length := total_length / float(segment_count)
	
	for i in range(segment_count):
		var y_offset := -segment_length * (i + 0.5)
		var xform := Transform3D(Basis.IDENTITY, Vector3(0, y_offset, 0))
		multimesh.set_instance_transform(i, xform)
