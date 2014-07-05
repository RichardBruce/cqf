# Third party library un-tar and compilation Makefile
CWD 	= $(shell pwd)

# SDL library for video output
# Might require libxv-dev installing
SDL2-$(SDL_VER) : SDL2-$(SDL_VER).tar.gz
	@echo "Unzipping SDL....."
	@tar xfz SDL2-$(SDL_VER).tar.gz
	@echo "Building SDL....."
	@(cd SDL2-$(SDL_VER)/ ; ./configure --silent --prefix $(CWD)/SDL2-$(SDL_VER) ; make -s ; make -s install)

# SDL TTF library for text rendering
SDL2_ttf-$(SDLTTF_VER) : SDL2_ttf-$(SDLTTF_VER).tar.gz freefont-ttf-$(FREEFONT_VER) freetype-$(FREETYPE_VER)
	@echo "Unzipping SDL TTF....."
	@tar xfz SDL2_ttf-$(SDLTTF_VER).tar.gz
	@echo "Building SDL TTF....."
	@(cd SDL2_ttf-$(SDLTTF_VER)/ ; ./configure --silent --prefix $(CWD)/SDL2_ttf-$(SDLTTF_VER) --with-freetype-prefix=$(CWD)/freetype-$(FREETYPE_VER) --with-sdl-prefix=$(CWD)/SDL2-$(SDL_VER); make -s ; make -s install)

# SDL image library for loading images to SDL
SDL2_image-$(SDLIMAGE_VER) : SDL2_image-$(SDLIMAGE_VER).tar.gz SDL2-$(SDL_VER)
	@echo "Unzipping SDL Image....."
	@tar xfz SDL2_image-$(SDLIMAGE_VER).tar.gz
	@echo "Building SDL Image....."
	@(cd SDL2_image-$(SDLIMAGE_VER)/ ; ./configure --silent --prefix $(CWD)/SDL2_image-$(SDLIMAGE_VER); make -s ; make -s install)

# FreeFont library for text rendering
freefont-ttf-$(FREEFONT_VER) : freefont-ttf-$(FREEFONT_VER).tar.gz
	@echo "Unzipping FreeFont....."
	@tar xfz freefont-ttf-$(FREEFONT_VER).tar.gz

# FreeType library for text rendering
freetype-$(FREETYPE_VER) : freetype-$(FREETYPE_VER).tar.gz
	@echo "Unzipping FreeType....."
	@tar xfz freetype-$(FREETYPE_VER).tar.gz
	@echo "Building FreeType....."
	@(cd freetype-$(FREETYPE_VER)/ ; ./configure --silent --prefix $(CWD)/freetype-$(FREETYPE_VER) ; make -s ; make -s install)

# TBB library for threading
tbb$(TBB_VER) : tbb$(TBB_VER)_src.tgz
	@echo "Unzipping TBB....."
	@tar xfz tbb$(TBB_VER)_src.tgz
	@echo "Building TBB..."
	@(cd tbb$(TBB_VER)/ ; make -s work_dir=./build/build )

# FFTW for optimised FFT implementation
fftw-$(FFTW_VER) : fftw-$(FFTW_VER).tar.gz
	@echo "Unzipping fftw-$(FFTW_VER)....."
	@tar xfz fftw-$(FFTW_VER).tar.gz
	@echo "Building fftw-$(FFTW_VER)....."
	@(cd fftw-$(FFTW_VER)/ ; ./configure --silent --prefix $(CWD)/fftw-$(FFTW_VER) --enable-float ; make -s ; make -s install)

# FFmpeg for optimised video encoding
ffmpeg : ffmpeg.tar.gz
	@echo "Unzipping ffmpeg....."
	@tar xfz ffmpeg.tar.gz
	@echo "Building ffmpeg....."
	@(cd ffmpeg/ ; ./configure --prefix=$(CWD)/ffmpeg --disable-indev=v4l --disable-indev=v4l2 --enable-memalign-hack ; make -s ; make -s install)
        
# SPD for test scenes
spd$(SPD_VER) : spd$(SPD_VER).tar.gz
	@echo "Unzipping spd$(SPD_VER)....."
	@tar xfz spd$(SPD_VER).tar.gz
	@echo "Building spd$(SPD_VER)....."
	@( cd spd$(SPD_VER) ; make -s all )

# Boost c++ libraries for c++ utilities
boost_$(BOOST_VER) : boost_$(BOOST_VER).tar.bz2
	@echo "Unzipping boost....."
	@tar xfj boost_$(BOOST_VER).tar.bz2
	@echo "Building boost....."
	@(cd boost_$(BOOST_VER)/tools/build/v2 ; ./bootstrap.sh ; ./b2 install --prefix=./ -d0)
	@(cd boost_$(BOOST_VER) ; ./bootstrap.sh variant=release --with-libraries=log,serialization,system,test --prefix=./ ; ./b2 -d0 -j8 install)
        
all :	SDL2-$(SDL_VER) SDL2_ttf-$(SDLTTF_VER) tbb$(TBB_VER) fftw-$(FFTW_VER) ffmpeg spd$(SPD_VER) boost_$(BOOST_VER)

# Clean targets
clean ::
	$(RM) -r SDL2-$(SDL_VER)
	$(RM) -r SDL2_ttf-$(SDLTTF_VER)
	$(RM) -r SDL2_image-$(SDLIMAGE_VER)
	$(RM) -r freefont-$(FREEFONT_VER)
	$(RM) -r freetype-$(FREETYPE_VER)
	$(RM) -r tbb$(TBB_VER)
	$(RM) -r fftw-$(FFTW_VER)
	$(RM) -r ffmpeg
	$(RM) -r spd$(SPD_VER)
	$(RM) -r boost_$(BOOST_VER)

spotless ::
	$(RM) -r SDL2-$(SDL_VER)
	$(RM) -r SDL2_ttf-$(SDLTTF_VER)
	$(RM) -r SDL2_image-$(SDLIMAGE_VER)
	$(RM) -r freefont-$(FREEFONT_VER)
	$(RM) -r freetype-$(FREETYPE_VER)
	$(RM) -r tbb$(TBB_VER)
	$(RM) -r fftw-$(FFTW_VER)
	$(RM) -r ffmpeg
	$(RM) -r spd$(SPD_VER)
	$(RM) -r boost_$(BOOST_VER)
