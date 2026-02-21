@tool
class_name WiggleGizmos extends EditorNode3DGizmoPlugin

func _init():
	create_material("main", Color(1, 1, 0))
	create_material("lines", Color(1, 1, 0, 0.5))
	create_material("axis_x", Color(1, 0, 0))
	create_material("axis_y", Color(0, 1, 0))
	create_material("axis_z", Color(0, 0, 1))
	create_material("constraint", Color(0.5, 0.5, 0.5, 0.5))

func _get_gizmo_name():
	return "WiggleGizmos"

func _has_gizmo(node):
	return node is WiggleShapeSphere or node is WiggleShapePlane or node is WiggleShapeLine or node is WiggleShapeBox

func _redraw(gizmo: EditorNode3DGizmo):
	gizmo.clear()
	var node = gizmo.get_node_3d()

	if node is WiggleShapeSphere:
		_draw_sphere(gizmo, node as WiggleShapeSphere)
	elif node is WiggleShapePlane:
		_draw_plane(gizmo, node as WiggleShapePlane)
	elif node is WiggleShapeLine:
		_draw_line(gizmo, node as WiggleShapeLine)
	elif node is WiggleShapeBox:
		_draw_box(gizmo, node as WiggleShapeBox)

func _draw_sphere(gizmo: EditorNode3DGizmo, sphere: WiggleShapeSphere):
	var lines = PackedVector3Array()
	var radius = sphere.Radius
	var steps = 32
	
	# Draw 3 circles for the sphere
	for i in range(steps):
		var a = float(i) / steps * TAU
		var b = float(i + 1) / steps * TAU
		
		# XY circle
		lines.push_back(Vector3(cos(a) * radius, sin(a) * radius, 0))
		lines.push_back(Vector3(cos(b) * radius, sin(b) * radius, 0))
		
		# XZ circle
		lines.push_back(Vector3(cos(a) * radius, 0, sin(a) * radius))
		lines.push_back(Vector3(cos(b) * radius, 0, sin(b) * radius))
		
		# YZ circle
		lines.push_back(Vector3(0, cos(a) * radius, sin(a) * radius))
		lines.push_back(Vector3(0, cos(b) * radius, sin(b) * radius))
		
	gizmo.add_lines(lines, get_material("main", gizmo))

func _draw_plane(gizmo: EditorNode3DGizmo, _plane: WiggleShapePlane):
	var lines = PackedVector3Array()
	var size = 2.0
	
	# Draw a square on the XZ plane (assuming Y is up/normal)
	lines.push_back(Vector3(-size, 0, -size))
	lines.push_back(Vector3(size, 0, -size))
	
	lines.push_back(Vector3(size, 0, -size))
	lines.push_back(Vector3(size, 0, size))
	
	lines.push_back(Vector3(size, 0, size))
	lines.push_back(Vector3(-size, 0, size))
	
	lines.push_back(Vector3(-size, 0, size))
	lines.push_back(Vector3(-size, 0, -size))
	
	# Draw normal
	lines.push_back(Vector3(0, 0, 0))
	lines.push_back(Vector3(0, 1, 0))
	
	gizmo.add_lines(lines, get_material("main", gizmo))

func _draw_line(gizmo: EditorNode3DGizmo, line: WiggleShapeLine):
	var lines = PackedVector3Array()
	lines.push_back(line.PointA)
	lines.push_back(line.PointB)
	gizmo.add_lines(lines, get_material("main", gizmo))

func _draw_box(gizmo: EditorNode3DGizmo, box: WiggleShapeBox):
	var lines = PackedVector3Array()
	var s = box.Size / 2.0
	
	# Bottom
	lines.push_back(Vector3(-s.x, -s.y, -s.z))
	lines.push_back(Vector3(s.x, -s.y, -s.z))
	lines.push_back(Vector3(s.x, -s.y, -s.z))
	lines.push_back(Vector3(s.x, -s.y, s.z))
	lines.push_back(Vector3(s.x, -s.y, s.z))
	lines.push_back(Vector3(-s.x, -s.y, s.z))
	lines.push_back(Vector3(-s.x, -s.y, s.z))
	lines.push_back(Vector3(-s.x, -s.y, -s.z))
	
	# Top
	lines.push_back(Vector3(-s.x, s.y, -s.z))
	lines.push_back(Vector3(s.x, s.y, -s.z))
	lines.push_back(Vector3(s.x, s.y, -s.z))
	lines.push_back(Vector3(s.x, s.y, s.z))
	lines.push_back(Vector3(s.x, s.y, s.z))
	lines.push_back(Vector3(-s.x, s.y, s.z))
	lines.push_back(Vector3(-s.x, s.y, s.z))
	lines.push_back(Vector3(-s.x, s.y, -s.z))
	
	# Verticals
	lines.push_back(Vector3(-s.x, -s.y, -s.z))
	lines.push_back(Vector3(-s.x, s.y, -s.z))
	lines.push_back(Vector3(s.x, -s.y, -s.z))
	lines.push_back(Vector3(s.x, s.y, -s.z))
	lines.push_back(Vector3(s.x, -s.y, s.z))
	lines.push_back(Vector3(s.x, s.y, s.z))
	lines.push_back(Vector3(-s.x, -s.y, s.z))
	lines.push_back(Vector3(-s.x, s.y, s.z))
	
	gizmo.add_lines(lines, get_material("main", gizmo))
