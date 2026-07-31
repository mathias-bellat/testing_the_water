FROM rocker/r-ver:4.4.0

RUN apt-get update && apt-get install -y \
    # Dépendances de base
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    # fs (libuv)
    libuv1-dev \
    # sf, terra
    libgdal-dev \
    gdal-bin \
    libgeos-dev \
    # units
    libudunits2-dev \
    # s2
    libabsl-dev \
    # exactextractr
    libgeos-dev \
    # ncdf4
    libnetcdf-dev \
    # igraph
    libglpk-dev \
    # systemfonts, textshaping
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    # png, svglite
    libpng-dev \
    # clipr
    libx11-dev \
    # SDMtune, dismo (Java)
    default-jdk \
    # rmarkdown, knitr
    pandoc \
    # git (remotes)
    git \
    # cmake (fs)
    cmake \
    && R CMD javareconf \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages('renv', repos='https://cloud.r-project.org')"

WORKDIR /app

# renv en premier pour profiter du cache Docker
COPY renv.lock renv.lock
COPY renv/ renv/
RUN R -e "renv::restore()"

COPY . .

CMD ["Rscript", "script.R"]
