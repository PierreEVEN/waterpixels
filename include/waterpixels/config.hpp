#pragma once

// Enable profiler recording and logging
#define ENABLE_PROFILER true

// Prompt intermediate images during generation
#define OUTPUT_DEBUG true

// The epsilon used in min value search (if abs(value - minValue) < E { take() })
#define WP_MARKER_EPSILON 0.0001

// Select the most central connected component as a source
#define PREFER_CELL_CENTER true

// Choose the Linf distance instead of L2 for spatial regularization
#define USE_LINF_REG_DISTANCE false
