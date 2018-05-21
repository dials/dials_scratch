import libtbx.load_env
import os
import platform
Import("env_etc")

env_etc.dials_scratch_dist = libtbx.env.dist_path("dials_scratch")
env_etc.dials_scratch_include = os.path.dirname(env_etc.dials_scratch_dist)
if (not env_etc.no_boost_python and hasattr(env_etc, "boost_adaptbx_include")):
    Import("env_no_includes_boost_python_ext")
    env = env_no_includes_boost_python_ext.Clone()
    env_etc.enable_more_warnings(env=env)
    env_etc.include_registry.append(
        env=env,
        paths=[
            env_etc.libtbx_include,
            env_etc.scitbx_include,
            env_etc.cctbx_include,
            env_etc.rstbx_include,
            env_etc.boost_include,
            env_etc.boost_adaptbx_include,
            env_etc.python_include,
            env_etc.dxtbx_include,
            env_etc.dials_include])
    env.SharedLibrary(
      target="#lib/dials_scratch_ext",
      source=["ext.cpp"])
	    
    env.Append(
                LIBS=env_etc.libm + [
                "scitbx_boost_python",
                "boost_python",
                "cctbx"])
    env.SConscript('idy/algorithms/spot_prediction/SConscript', exports={ 'env' : env })
    env.SConscript('dgw/cctbx_cpp_examples/SConscript', exports={ 'env' : env })
    env.SConscript('jmp/SConscript', exports={ 'env' : env })
