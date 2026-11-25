R"csjson({
  "profile": {
    "property_type": "EnumProperty",
    "name": "Profile",
    "items": [
      [
        "auto",
        "Auto",
        "Choose from geometric panels, revolved parts and curved surfaces"
      ],
      [
        "hard_surface",
        "Hard surface",
        "Seams on hard edges, straight sides where the topology is a grid"
      ],
      [
        "organic",
        "Organic",
        "No feature loops, smoothed contours that keep the curves"
      ],
      [
        "custom",
        "Custom",
        "Whatever is set below"
      ]
    ],
    "update": {
      "expression": "_profile_update"
    }
  },
  "generate": {
    "property_type": "BoolProperty",
    "name": "Generate seams",
    "description": "Trace the quad layout and merge patches under the distortion budget. Off: only refine the existing seams and clean the contours"
  },
  "protect_existing": {
    "property_type": "BoolProperty",
    "name": "Protect existing seams",
    "description": "Seams already on the mesh are kept as locked cuts"
  },
  "feature_angle": {
    "property_type": "FloatProperty",
    "name": "Feature angle",
    "min": 0.0,
    "max": 180.0,
    "description": "Edges sharper than this (degrees) become seams and are never merged away; 0 disables"
  },
  "refine_paths": {
    "property_type": "BoolProperty",
    "name": "Refine seam paths",
    "description": "Re-route open seam chains inside a corridor so they stop hopping between quad rows; anchors, junctions and locked seams are kept"
  },
  "contour_mode": {
    "property_type": "EnumProperty",
    "name": "Contours",
    "items": [
      [
        "auto",
        "Auto",
        "Straighten arcs that are edge loops and already nearly straight, smooth the rest"
      ],
      [
        "panels",
        "Panels",
        "Straighten every arc between anchors"
      ],
      [
        "organic",
        "Organic",
        "Smooth short oscillations, keep the overall curve"
      ],
      [
        "off",
        "Off",
        "Leave the unwrap's contours as they are"
      ]
    ]
  },
  "budget": {
    "property_type": "EnumProperty",
    "name": "Distortion budget",
    "items": [
      [
        "tight",
        "Tight",
        "P95 anisotropy 1.3, area 1.4"
      ],
      [
        "normal",
        "Normal",
        "P95 anisotropy 1.6, area 1.8"
      ],
      [
        "loose",
        "Loose",
        "P95 anisotropy 2.0, area 2.4"
      ]
    ]
  },
  "min_island": {
    "property_type": "FloatProperty",
    "name": "Absorb islands below %",
    "min": 0.0,
    "max": 10.0,
    "description": "Islands smaller than this share of the surface are merged into a neighbour under a relaxed budget"
  },
  "solver_iterations": {
    "property_type": "IntProperty",
    "name": "Solver iterations",
    "min": 1,
    "max": 50,
    "description": "Minimum-stretch iterations of Blender's unwrap"
  },
  "min_island_faces": {
    "property_type": "IntProperty",
    "name": "Minimum island faces",
    "min": 0,
    "description": "Absorb smaller islands into a validated neighbour, except whole components and boundaries protected by the artist"
  },
  "target_enabled": {
    "property_type": "BoolProperty",
    "name": "Target island count",
    "description": "Search up to five distortion budgets; an unreachable target is reported"
  },
  "target_min": {
    "property_type": "IntProperty",
    "name": "Minimum",
    "min": 1
  },
  "target_max": {
    "property_type": "IntProperty",
    "name": "Maximum",
    "min": 1
  },
  "deadline_s": {
    "property_type": "FloatProperty",
    "name": "Time budget (seconds)",
    "min": 0.0,
    "description": "Return a validated partial layout when time runs out; 0 is unlimited. An active solve finishes before the budget is checked again"
  },
  "hide_from": {
    "property_type": "EnumProperty",
    "name": "Hide seams from",
    "items": [
      [
        "none",
        "All views equally",
        "Use overall surface visibility"
      ],
      [
        "view",
        "Current view",
        "Prefer seams away from the current 3D view"
      ],
      [
        "custom",
        "Custom direction",
        "Direction toward the viewer in object coordinates"
      ],
      [
        "+X",
        "+X object axis",
        "Keep seams away from the +X side"
      ],
      [
        "-X",
        "-X object axis",
        "Keep seams away from the -X side"
      ],
      [
        "+Y",
        "+Y object axis",
        "Keep seams away from the +Y side"
      ],
      [
        "-Y",
        "-Y object axis",
        "Keep seams away from the -Y side"
      ],
      [
        "+Z",
        "+Z object axis",
        "Keep seams away from the +Z side"
      ],
      [
        "-Z",
        "-Z object axis",
        "Keep seams away from the -Z side"
      ]
    ]
  },
  "hide_direction": {
    "property_type": "FloatVectorProperty",
    "name": "Viewing direction",
    "size": 3,
    "subtype": "DIRECTION",
    "description": "Direction from the object toward the viewer, in object coordinates"
  },
  "show_advanced": {
    "property_type": "BoolProperty",
    "name": "Advanced"
  },
  "w_turn": {
    "property_type": "FloatProperty",
    "name": "Turn weight",
    "min": 0.0,
    "max": 20.0
  },
  "w_flow": {
    "property_type": "FloatProperty",
    "name": "Flow weight",
    "min": 0.0,
    "max": 20.0
  },
  "w_drift": {
    "property_type": "FloatProperty",
    "name": "Drift weight",
    "min": 0.0,
    "max": 20.0
  },
  "corridor_rings": {
    "property_type": "IntProperty",
    "name": "Corridor rings",
    "min": 1,
    "max": 6
  },
  "smooth_radius": {
    "property_type": "FloatProperty",
    "name": "Smooth radius",
    "min": 0.5,
    "max": 12.0,
    "description": "In local edge lengths"
  },
  "max_displacement": {
    "property_type": "FloatProperty",
    "name": "Max displacement",
    "min": 0.1,
    "max": 6.0,
    "description": "Band for smoothed contours, in edge lengths"
  },
  "min_fusion": {
    "property_type": "FloatProperty",
    "name": "Min fusion",
    "min": 0.0,
    "max": 0.9,
    "description": "A merge must remove at least this share of the smaller island's outline"
  },
  "pack_after": {
    "property_type": "EnumProperty",
    "name": "Then pack with",
    "items": [
      [
        "none",
        "Nothing",
        "Leave the islands where the unwrap put them"
      ],
      [
        "blender",
        "Blender",
        "bpy.ops.uv.pack_islands"
      ]
    ]
  },
  "keep_orientation": {
    "property_type": "BoolProperty",
    "name": "Keep island orientation"
  },
  "pack_width": {
    "property_type": "IntProperty",
    "name": "Texture width",
    "min": 1
  },
  "pack_height": {
    "property_type": "IntProperty",
    "name": "Texture height",
    "min": 1
  },
  "pack_padding": {
    "property_type": "FloatProperty",
    "name": "Island gap (texels)",
    "min": 0
  },
  "pack_border": {
    "property_type": "FloatProperty",
    "name": "Tile border (texels)",
    "min": 0
  },
  "pack_tiles": {
    "property_type": "StringProperty",
    "name": "UDIM tiles",
    "description": "Comma-separated UDIM tiles, for example 1001,1002"
  },
  "last_status": {
    "property_type": "StringProperty",
    "name": "Last status"
  },
  "last_report": {
    "property_type": "StringProperty",
    "name": "Last report"
  },
  "policy": {
    "property_type": "EnumProperty",
    "name": "Seam policy",
    "description": "Native baseline, the default, or the enhanced proposal pass after it",
    "items": [
      [
        "native_baseline",
        "Native baseline",
        "The validated baseline layout: every correction of the native engine, no proposal pass"
      ],
      [
        "enhanced_v1",
        "Enhanced",
        "The native baseline first, then bounded, validated proposals (sliding, relief, regions, rings, cavities, shells, symmetry, spur removal); the baseline is kept unless a proposal ranks strictly better"
      ]
    ],
    "group": "enhanced"
  },
  "show_enhanced": {
    "property_type": "BoolProperty",
    "name": "Enhanced settings",
    "group": "enhanced"
  },
  "enhanced_sliding": {
    "property_type": "BoolProperty",
    "name": "Slide tube slits",
    "description": "Propose moving the slit of a tube-like part to its least exposed side along the rims; applied only when the validated result ranks better",
    "group": "enhanced"
  },
  "max_slide_distance": {
    "property_type": "FloatProperty",
    "name": "Max slide distance",
    "description": "Object-space distance a slit end may move along a rim; 0 is unlimited",
    "min": 0.0,
    "group": "enhanced"
  },
  "enhanced_relief": {
    "property_type": "BoolProperty",
    "name": "Relief slits",
    "description": "Short slits at distortion hotspots of over-budget islands, before any broad split",
    "group": "enhanced"
  },
  "enhanced_regions": {
    "property_type": "BoolProperty",
    "name": "Preference-aware regions",
    "description": "Re-partition over-budget islands along welcome seams (painted allow, creases), away from avoided ones",
    "group": "enhanced"
  },
  "enhanced_rings": {
    "property_type": "BoolProperty",
    "name": "Constriction rings",
    "description": "Closed cuts around necks and creases, with their companion cuts",
    "group": "enhanced"
  },
  "ring_confidence": {
    "property_type": "IntProperty",
    "name": "Ring confidence %",
    "description": "Rings with less constriction or crease evidence are reported, not proposed (percent of the full evidence)",
    "min": 0,
    "max": 100,
    "group": "enhanced"
  },
  "enhanced_cavities": {
    "property_type": "BoolProperty",
    "name": "Cavity rings",
    "description": "Rings and interior slits around occluded recesses",
    "group": "enhanced"
  },
  "cavity_confidence": {
    "property_type": "IntProperty",
    "name": "Cavity confidence %",
    "description": "Recesses with less occlusion evidence are reported, not proposed (percent of the full evidence)",
    "min": 0,
    "max": 100,
    "group": "enhanced"
  },
  "enhanced_shells": {
    "property_type": "BoolProperty",
    "name": "Thin-shell rims",
    "description": "Front/back splits of thin parts along their rims",
    "group": "enhanced"
  },
  "shell_confidence": {
    "property_type": "IntProperty",
    "name": "Shell confidence %",
    "description": "Thin parts with less front/back evidence are reported, not proposed (percent of the full evidence)",
    "min": 0,
    "max": 100,
    "group": "enhanced"
  },
  "enhanced_spur_removal": {
    "property_type": "BoolProperty",
    "name": "Remove unnecessary spurs",
    "description": "Remove seam dead ends whose removal passes validation; necessary relief ends stay",
    "group": "enhanced"
  },
  "symmetry_mode": {
    "property_type": "EnumProperty",
    "name": "Symmetry",
    "description": "Whether seam moves are mirrored",
    "items": [
      [
        "asymmetric",
        "Asymmetric",
        "View-optimal: moves are not paired"
      ],
      [
        "strict",
        "Strict",
        "Every move is paired with its mirror image in one validated operation; needs an exact mirror correspondence"
      ],
      [
        "approximate",
        "Approximate (report only)",
        "Symmetry metrics are reported; no paired moves and no hard edge equality"
      ]
    ],
    "group": "enhanced"
  },
  "symmetry_on_conflict": {
    "property_type": "EnumProperty",
    "name": "On a mirror conflict",
    "description": "What strict symmetry does when the mirror half breaks a lock, pin, anchor or painted constraint",
    "items": [
      [
        "reject_pair",
        "Reject the pair",
        "Keep the current layout"
      ],
      [
        "disable_for_proposal",
        "Evaluate unpaired",
        "Evaluate the move alone, labelled as unpaired in the report"
      ]
    ],
    "group": "enhanced"
  },
  "symmetry_plane": {
    "property_type": "EnumProperty",
    "name": "Mirror plane",
    "description": "Detected, or an object-axis plane at the offset",
    "items": [
      [
        "detect",
        "Detect",
        "Find the mirror plane of the mesh"
      ],
      [
        "X",
        "X",
        "The plane x = offset, in object coordinates"
      ],
      [
        "Y",
        "Y",
        "The plane y = offset, in object coordinates"
      ],
      [
        "Z",
        "Z",
        "The plane z = offset, in object coordinates"
      ]
    ],
    "group": "enhanced"
  },
  "symmetry_offset": {
    "property_type": "FloatProperty",
    "name": "Plane offset",
    "description": "Position of an axis mirror plane along its axis, in object coordinates",
    "group": "enhanced"
  },
  "scale_mode": {
    "property_type": "EnumProperty",
    "name": "Distances",
    "description": "How the enhanced features measure distances and sizes",
    "items": [
      [
        "physical",
        "Physical",
        "Object-space lengths and areas, frozen at the start of the pass"
      ],
      [
        "hops",
        "Edge hops (legacy)",
        "Hop rings, face counts and fixed step counts, as the native baseline measures"
      ]
    ],
    "group": "enhanced"
  },
  "reference_length": {
    "property_type": "FloatProperty",
    "name": "Reference length",
    "description": "Object-space length unit of the physical distances; 0 uses the median edge length. Give the same value to compare different tessellations",
    "min": 0.0,
    "group": "enhanced"
  },
  "required_attribute": {
    "property_type": "StringProperty",
    "name": "Required seams",
    "description": "EDGE BOOLEAN attribute of edges that must be seams; empty disables",
    "group": "enhanced"
  },
  "forbidden_attribute": {
    "property_type": "StringProperty",
    "name": "Forbidden seams",
    "description": "EDGE BOOLEAN attribute of edges that must never become seams; empty disables",
    "group": "enhanced"
  },
  "allow_attribute": {
    "property_type": "StringProperty",
    "name": "Welcome seams",
    "description": "POINT FLOAT attribute in [0, 1] where seams are welcome; empty disables",
    "group": "enhanced"
  },
  "avoid_attribute": {
    "property_type": "StringProperty",
    "name": "Avoided seams",
    "description": "POINT FLOAT attribute in [0, 1] where seams should be avoided; empty disables",
    "group": "enhanced"
  },
  "max_candidates_per_chart": {
    "property_type": "IntProperty",
    "name": "Candidates per island",
    "description": "Complete alternative layouts evaluated per island in one pass",
    "min": 1,
    "max": 1024,
    "group": "enhanced"
  },
  "max_solves_per_chart": {
    "property_type": "IntProperty",
    "name": "Solves per island",
    "description": "Unwraps of proposals per island in one pass",
    "min": 1,
    "max": 4096,
    "group": "enhanced"
  },
  "prune_dominated": {
    "property_type": "BoolProperty",
    "name": "Skip provably worse",
    "description": "Refuse, before unwrapping, a proposal whose seams alone prove it cannot beat the kept layout",
    "group": "enhanced"
  },
  "last_messages": {
    "property_type": "StringProperty",
    "name": "Last messages",
    "group": "enhanced"
  }
})csjson"
