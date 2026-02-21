@tool
class_name ChainWiggler extends MultiMeshInstance3D

const COMPLIANCE_FACTOR : float = 0.00001

@export var link_count: int = 5:
	set(value):
		link_count = max(1, value)

@export var total_length: float = 2.0
@export var link_width: float = 0.1

@export var link_mass: float = 0.1:
	set(value):
		link_mass = value
		for p_id in _link_particles:
			if p_id != -1: WiggleKit.particle_set_mass(p_id, value)

@export var link_damping: float = 0.5:
	set(value):
		link_damping = value
		for p_id in _link_particles:
			if p_id != -1: WiggleKit.particle_set_damping(p_id, value)
		for p_id in _aux_particles:
			if p_id != -1: WiggleKit.particle_set_damping(p_id, value)

@export var gravity: Vector3 = Vector3(0, -10, 0):
	set(value):
		gravity = value
		for p_id in _link_particles:
			if p_id != -1: WiggleKit.particle_set_gravity(p_id, value)
		for p_id in _aux_particles:
			if p_id != -1: WiggleKit.particle_set_gravity(p_id, value)

@export var stretch_compliance: float = 0.0
@export var twist_compliance: float = 1.0

@export var Colliders: Array[WiggleShapeBase] = []
@export var ShapeDistanceConstraints: Array[WiggleShapeBase] = []

@export var end_attachment: Node3D = null

# Fixed anchor particle (mass 0)
var _anchor_particle: int = -1
var _anchor_aux_particles: Array[int] = []  # 3 auxiliary particles for anchor triangle

# Chain link particles (main axis)
var _link_particles: Array[int] = []
# Auxiliary particles for each link (3 per link forming a triangle)
var _aux_particles: Array[int] = []

var _distance_constraints: Array[int] = []
var _tetrahedral_constraints: Array[int] = []


func _enter_tree() -> void:
	if Engine.is_editor_hint():
		_setup_multimesh()
		return
	
	_setup_multimesh()
	_create_particles()
	_create_constraints()
	_register_colliders()
	await get_tree().process_frame
	if end_attachment:
		WiggleKit.particle_set_position(_link_particles[link_count - 1], end_attachment.global_position)
		var all_particles : Array
		all_particles.append_array(_link_particles)
		all_particles.append_array(_aux_particles)
		WiggleKit.warmup(200, all_particles)


func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if WiggleKit:
		for c_id in _tetrahedral_constraints:
			if c_id != -1: WiggleKit.tetrahedral_constraint_free(c_id)
		_tetrahedral_constraints.clear()
		
		for c_id in _distance_constraints:
			if c_id != -1: WiggleKit.distance_constraint_free(c_id)
		_distance_constraints.clear()
		
		if _anchor_particle != -1: WiggleKit.particle_free(_anchor_particle)
		for p_id in _anchor_aux_particles:
			if p_id != -1: WiggleKit.particle_free(p_id)
		_anchor_aux_particles.clear()
		
		for p_id in _link_particles:
			if p_id != -1: WiggleKit.particle_free(p_id)
		_link_particles.clear()
		
		for p_id in _aux_particles:
			if p_id != -1: WiggleKit.particle_free(p_id)
		_aux_particles.clear()


func _setup_multimesh() -> void:
	if multimesh == null:
		multimesh = MultiMesh.new()
		multimesh.transform_format = MultiMesh.TRANSFORM_3D
	multimesh.instance_count = link_count


func _create_particles() -> void:
	var link_length := total_length / float(link_count)
	var half_width := link_width / 2.0
	
	# Create anchor (fixed) particle at origin
	_anchor_particle = WiggleKit.particle_create(global_position, 0, Vector3.ZERO, 0)
	
	# Create 3 auxiliary particles for anchor forming a triangle in XZ plane
	var anchor_aux_offsets := [
		Vector3(half_width, 0, 0),
		Vector3(-half_width * 0.5, 0, half_width * 0.866),
		Vector3(-half_width * 0.5, 0, -half_width * 0.866)
	]
	for offset in anchor_aux_offsets:
		var p_id = WiggleKit.particle_create(to_global(offset), 0, Vector3.ZERO, 0)
		_anchor_aux_particles.append(p_id)
	
	# Create chain link particles
	for i in range(link_count):
		var y_offset := -link_length * (i + 1)
		var pos := to_global(Vector3(0, y_offset, 0))
		var mass := link_mass
		if end_attachment and i == link_count - 1: mass = 0.0
		var p_id = WiggleKit.particle_create(pos, mass, gravity, link_damping)
		_link_particles.append(p_id)
		
		# Create 3 auxiliary particles per link (triangle)
		var aux_offsets := [
			Vector3(half_width, y_offset, 0),
			Vector3(-half_width * 0.5, y_offset, half_width * 0.866),
			Vector3(-half_width * 0.5, y_offset, -half_width * 0.866)
		]
		for offset in aux_offsets:
			var aux_pos := to_global(offset)
			var aux_id = WiggleKit.particle_create(aux_pos, link_mass * 0.1, gravity, link_damping)
			_aux_particles.append(aux_id)


func _create_constraints() -> void:
	var stretch_comp := stretch_compliance * COMPLIANCE_FACTOR
	var twist_comp := twist_compliance * COMPLIANCE_FACTOR
	
	# Anchor triangle internal constraints (rigid)
	for i in range(3):
		var next := (i + 1) % 3
		_distance_constraints.append(WiggleKit.distance_constraint_create(
			_anchor_aux_particles[i], _anchor_aux_particles[next], 0))
		_distance_constraints.append(WiggleKit.distance_constraint_create(
			_anchor_particle, _anchor_aux_particles[i], 0))
	
	# For each link
	for i in range(link_count):
		var main_p := _link_particles[i]
		var aux_base := i * 3
		
		# Internal triangle constraints (keep link rigid)
		for j in range(3):
			var next := (j + 1) % 3
			_distance_constraints.append(WiggleKit.distance_constraint_create(
				_aux_particles[aux_base + j], _aux_particles[aux_base + next], 0))
			_distance_constraints.append(WiggleKit.distance_constraint_create(
				main_p, _aux_particles[aux_base + j], 0))
		
		# Connect to previous link (or anchor)
		var prev_main: int
		var prev_aux_base: int
		var prev_aux: Array[int] = []
		
		if i == 0:
			prev_main = _anchor_particle
			prev_aux = []
			prev_aux.append(_anchor_aux_particles[0])
			prev_aux.append(_anchor_aux_particles[1])
			prev_aux.append(_anchor_aux_particles[2])
		else:
			prev_main = _link_particles[i - 1]
			var pb := (i - 1) * 3
			prev_aux = []
			prev_aux.append(_aux_particles[pb])
			prev_aux.append(_aux_particles[pb + 1])
			prev_aux.append(_aux_particles[pb + 2])
		
		# Distance constraint between main particles (chain backbone)
		_distance_constraints.append(WiggleKit.distance_constraint_create(
			prev_main, main_p, stretch_comp))
		
		# Cross-link distance constraints for stability
		for j in range(3):
			_distance_constraints.append(WiggleKit.distance_constraint_create(
				prev_aux[j], _aux_particles[aux_base + j], stretch_comp))
		
		# Tetrahedral constraints for twist resistance
		# Create tetrahedra spanning both links
		for j in range(3):
			var next := (j + 1) % 3
			_tetrahedral_constraints.append(WiggleKit.tetrahedral_constraint_create(
				prev_main, prev_aux[j], _aux_particles[aux_base + j], _aux_particles[aux_base + next], twist_comp))


func _register_colliders() -> void:
	for p_id in _link_particles:
		for collider in Colliders:
			if collider: collider.add_particle_to_collider(p_id)
		for sdf in ShapeDistanceConstraints:
			if sdf: sdf.add_particle_to_sdf_constraint(p_id)
	
	for p_id in _aux_particles:
		for collider in Colliders:
			if collider: collider.add_particle_to_collider(p_id)
		for sdf in ShapeDistanceConstraints:
			if sdf: sdf.add_particle_to_sdf_constraint(p_id)


func _physics_process(_delta: float) -> void:
	if Engine.is_editor_hint():
		_update_editor_preview()
		return
	
	if _link_particles.size() < link_count:
		return
	
	# Update anchor position to follow node
	WiggleKit.particle_set_position(_anchor_particle, global_position)
	
	# Update anchor auxiliary particles
	var half_width := link_width / 2.0
	var anchor_aux_offsets := [
		Vector3(half_width, 0, 0),
		Vector3(-half_width * 0.5, 0, half_width * 0.866),
		Vector3(-half_width * 0.5, 0, -half_width * 0.866)
	]
	for i in range(3):
		WiggleKit.particle_set_position(_anchor_aux_particles[i], to_global(anchor_aux_offsets[i]))
	
	# Update end attachment if set
	if end_attachment and is_instance_valid(end_attachment):
		var last_idx := link_count - 1
		WiggleKit.particle_set_position(_link_particles[last_idx], end_attachment.global_position)
	
	# Update multimesh instances
	_update_multimesh()


func _update_multimesh() -> void:
	if multimesh == null or _link_particles.size() < link_count:
		return
	
	var link_length := total_length / float(link_count)
	
	for i in range(link_count):
		var main_pos := WiggleKit.particle_get_position(_link_particles[i])
		var aux_base := i * 3
		var aux0 := WiggleKit.particle_get_position(_aux_particles[aux_base])
		var aux1 := WiggleKit.particle_get_position(_aux_particles[aux_base + 1])
		var aux2 := WiggleKit.particle_get_position(_aux_particles[aux_base + 2])
		
		# Get previous link position for direction
		var prev_pos: Vector3
		if i == 0:
			prev_pos = WiggleKit.particle_get_position(_anchor_particle)
		else:
			prev_pos = WiggleKit.particle_get_position(_link_particles[i - 1])
		
		# Calculate basis from particle positions
		var y_axis := (prev_pos - main_pos).normalized()
		var x_axis := (aux0 - main_pos).normalized()
		
		if y_axis.length_squared() < 0.0001:
			y_axis = Vector3.UP
		if x_axis.length_squared() < 0.0001:
			x_axis = Vector3.RIGHT
		
		# Orthonormalize
		x_axis = (x_axis - y_axis * x_axis.dot(y_axis)).normalized()
		var z_axis := x_axis.cross(y_axis).normalized()
		
		var basis := Basis(x_axis, y_axis, z_axis).orthonormalized()
		
		# Position at midpoint between this and previous link
		var mid_pos := (prev_pos + main_pos) / 2.0
		
		var xform := Transform3D(basis, mid_pos)
		multimesh.set_instance_transform(i, global_transform.affine_inverse() * xform)


func _update_editor_preview() -> void:
	if multimesh == null:
		return
	
	var link_length := total_length / float(link_count)
	
	for i in range(link_count):
		var y_offset := -link_length * (i + 0.5)
		var xform := Transform3D(Basis.IDENTITY, Vector3(0, y_offset, 0))
		multimesh.set_instance_transform(i, xform)
