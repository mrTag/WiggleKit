@tool
class_name WiggleShapeLine extends WiggleShapeBase

@export var PointA : Vector3 = Vector3(0, 0, 0):
	set(value):
		PointA = value
		update_gizmos()
@export var PointB : Vector3 = Vector3(0, 1, 0):
	set(value):
		PointB = value
		update_gizmos()

@export_group("SDF Constraint")
@export var SDF_Compliance : float = 0.0
@export var SDF_Distance : float = 0.0

var _sdf_id : int = -1

func update_gizmos():
	pass

func add_particle_to_collider(_particle_index:int):
	# Line doesn't support collider
	pass

func add_particle_to_sdf_constraint(particle_index:int):
	if Engine.is_editor_hint():
		return
	if _sdf_id == -1:
		_sdf_id = WiggleKit.line_sdf_create(to_global(PointA), to_global(PointB), SDF_Distance, SDF_Compliance, false)
	WiggleKit.line_sdf_add_particle(_sdf_id, particle_index)

func _process(_delta):
	if Engine.is_editor_hint():
		return
	if _sdf_id != -1:
		WiggleKit.line_sdf_set_position_a(_sdf_id, to_global(PointA))
		WiggleKit.line_sdf_set_position_b(_sdf_id, to_global(PointB))

func _exit_tree() -> void:
	if Engine.is_editor_hint():
		return
	if _sdf_id != -1:
		WiggleKit.line_sdf_free(_sdf_id)
		_sdf_id = -1
