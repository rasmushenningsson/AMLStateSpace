function download_samples(; dest_path=nothing)
	if dest_path === nothing
		dest_path = "data/samples"
		isdir(dest_path) || mkdir(dest_path)
	end

	figshare_path = Downloads.download("https://api.figshare.com/v2/articles/23715648")
	figshare_info = JSON.parsefile(figshare_path)
	figshare_files = figshare_info["files"]

	# name => md5
	sample_md5 = Pair{String,String}[
		"AML7.h5" => "5c2a19014832f96bba5881866a7f2119",
		"AML9.h5" => "f410d8daa9ba0a10548cd2efb8fb024b",
		"AML10.h5" => "676b645852b95fd221e56dc9b7cd1c49",
		"AML21.h5" => "8f6f1a9f3ad551ae29c65c9be1fd9f11",
		"AML23.h5" => "b0e4c17d7aa7143dc8a9dfe9eaec450a",
		"AML24.h5" => "9c59c6c9dd719020d6ca700c210757ae",
		"AML25.h5" => "21dbae564fad2e79708864ad5fce7874",
		"AML27.h5" => "11714c580e637c9f450a004002d31144",
		"AML28.h5" => "d408edc41e09685f57f103cbe6950e50",
		"AML32.h5" => "086f522ec749c25a086908ae6ff4078b",
		"AML33.h5" => "438dc76f959607e48681346e6051d1bc",
		"AML34.h5" => "17cb78d8b6c21bf5437f72e6b5c286e7",
		"AML37.h5" => "797416ef308d89d98ed9277ee32ce4a8",
		"AML48.h5" => "03a9da9750723bef5b4125eedb3cb983",
		"AML55.h5" => "12998f043f1d911c954726a32eae5840",
		"AML61.h5" => "55defeef5daa04da92f4877c83084045",
		"AML62.h5" => "b26585c7e4d89afbafbf65a769ae797d",
		"AML79.h5" => "7cdf532dfb509febb02d47faca0d05ab",
		"AML80.h5" => "903dd911ee776ccc0b2ca032a992a146",
		"AML83.h5" => "26c6327560af6a38740dda7a43e7ee7c",
		"AML85D.h5" => "3fcbe29ae225d2da0643453a8fa9b9d3",
		"AML85R.h5" => "e4f339626020c61e0b8295127df10b3e",
		"AML97.h5" => "e4aa5e7817bed9f09565fcc3e6721985",
		"AML104.h5" => "6654e42969fcc47194e05d641243976e",
		"AML105.h5" => "11234f9e7eb2182d9d3f4ba8e51a9d23",
		"AML110.h5" => "a7b08aa30c57ebefb16dc0d156827963",
		"AML111.h5" => "733e4c88167751598c344ad797da55c8",
		"AML117.h5" => "18c211f5f630d5f6c33bd82060759559",
		"AML123.h5" => "373e020cc1f6e7c75d6190dd7072f3c8",
		"AML124.h5" => "4d71d799cbecd7fea6d0280ef51bd294",
		"AML126.h5" => "42174579a0ba0fee0731c8c19651c23e",
		"AML136.h5" => "8809396668a7c929e5fbc958ecfc960e",
		"AML138.h5" => "416d7cff17e0e6a85d22267f83a0461e",
		"AML151.h5" => "8e9336f98eea4922a7b8f4be4f9935c6",
		"AML155.h5" => "84baacecdc209ed63898f9f2e61ef714",
		"AML157.h5" => "7a8e8e9e5b0d37a9007bd431bbb15c1e",
		"AML161.h5" => "e4ff45c00ff02a0049b34d23b5e0f2bd",
		"AML172.h5" => "97207962e179f0ba16487277ee8dd770",
		"NBM4-MNC.h5" => "6b47ef374d628310227ee4882309607f",
		"NBM7-MNC.h5" => "0722d94f14d2f62f79465f33dbd3363c",
		"NBM8-CD34.h5" => "226032d4e38f656e5f845914af9961bd",
		"NBM8-MNC.h5" => "e60be07df3828c9c675d5c1aa9e642f5",
		"NBM10-CD34.h5" => "1dba60200dbc4ad8fbaf23bdbe494118",
		"NBM10-MNC.h5" => "f9be34eea9a7776de70d96803402af15",
		"NBM11-CD34.h5" => "7ae94440bf3a1b1ada716a3e8006efe6",
		"NBM11-MNC.h5" => "f69b657c061d6d4b263b06e12e7f2412",
	]

	# name => url
	figshare_dict = Dict{String,String}()

	for file_dict in figshare_files
		name = get(file_dict, "name", nothing)
		name === nothing && continue
		download_url = get(file_dict, "download_url", nothing)
		download_url === nothing && continue
		figshare_dict[name] = download_url
	end


	for (name, expected_md5) in sample_md5
		file_path = joinpath(dest_path, name)
		if isfile(file_path)
			new_md5 = bytes2hex(open(MD5.md5,file_path))
			new_md5 === expected_md5 && continue # File exists and md5 is correct
			@warn "Mismatching md5 sum for sample $name, expected $expected_md5, but got $new_md5. Downloading again."
			rm(file_path)
		end

		# Download the file
		figshare_url = get(figshare_dict, name, nothing)
		@assert figshare_url !== nothing "Could not find sample $name in figshare."

		@info "Downloading $name"
		Downloads.download(figshare_url, file_path)

		new_md5 = bytes2hex(open(MD5.md5,file_path))
		@assert new_md5 === expected_md5 "Fatal error, downloaded file does not have the correct md5 sum, expected $expected_md5, but got $new_md5."
	end

	@info "All samples are downloaded and md5 sums verified"
end
