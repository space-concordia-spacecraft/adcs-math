#pragma once

#include "Constants.h"

#include "glm/glm.hpp"

using glm::vec3;
using glm::vec4;
using glm::mat3;
using glm::mat4;

#define CONST_DEG_TO_RAD CONST_PI / 180.0

typedef enum
{
    e80,
    e96,
    e00a,
    e00b,
    e00cio
} eOpt;

namespace adcs {
    bool isEclipse(vec3 se, vec3 pe);
    double gstime(double jdut1);
}