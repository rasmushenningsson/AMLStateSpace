getdelim(filepath::String) = lowercase(splitext(filepath)[2])==".csv" ? ',' : '\t'

set_samples_path(path::String) = @set_preferences!("samples_path"=>expanduser(path))
set_cache_path(path::String) = @set_preferences!("cache_path"=>expanduser(path))

get_samples_path() = @something @load_preference("samples_path") "data/samples"
get_annotations_path() = "data/annotations"
get_cache_path() = @something @load_preference("cache_path") ".cache"
