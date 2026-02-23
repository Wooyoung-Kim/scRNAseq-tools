FROM mambaorg/micromamba:1.5.8

COPY environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml -n scrnaseq && \
    micromamba clean -a -y

ENV PATH=/opt/conda/envs/scrnaseq/bin:$PATH
WORKDIR /workspace

COPY . /workspace
RUN pip install -e .

CMD ["scRNAseq-tools", "--help"]
