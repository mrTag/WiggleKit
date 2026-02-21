@tool
class_name WiggleShapePlane extends WiggleShapeBase

@export_group("Collider")
@export_range(0.0, 1.0, 0.05) var Bounciness : float = 1
@export_range(0.0, 1.0, 0.05) var Friction : float = 0

@export_group("SDF Constraint")
@export var SDF_Compliance : float = 0.0
@export var SDF_Distance : float = 0.0

var _collider_id : int = -1
var _sdf_id : int = -1

func add_particle_to_collider(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _collider_id == -1:
		if is_inside_tree():
			_collider_id = WiggleKit.plane_collider_create(global_position, global_basis.y, Bounciness, Friction)
		else:
			_collider_id = WiggleKit.plane_collider_create(position, basis.y, Bounciness, Friction)
	WiggleKit.plane_collider_add_particle(_collider_id, particle_index)

func add_particle_to_sdf_constraint(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _sdf_id == -1:
		_sdf_id = WiggleKit.plane_sdf_create(global_position, global_basis.y, SDF_Distance, SDF_Compliance)
	WiggleKit.plane_sdf_add_particle(_sdf_id, particle_index)

func _physics_process(_delta):
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.plane_collider_set_position(_collider_id, global_position)
		WiggleKit.plane_collider_set_normal(_collider_id, global_basis.y)
	if _sdf_id != -1:
		WiggleKit.plane_sdf_set_position(_sdf_id, global_position)
		WiggleKit.plane_sdf_set_normal(_sdf_id, global_basis.y)

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if _collider_id != -1:
		WiggleKit.plane_collider_free(_collider_id)
		_collider_id = -1
	if _sdf_id != -1:
		WiggleKit.plane_sdf_free(_sdf_id)
		_sdf_id = -1
