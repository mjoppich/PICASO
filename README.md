# PICASO


## Creating environment

    conda create --name regnetworks -c conda-forge python=3.11


    conda activate regnetworks

    conda install -c conda-forge python-igraph ipykernel
    python -m ipykernel install --user --name=regnetworks

    pip install scanpy networkx scikit-learn natsort seaborn matplotlib leidenalg markov_clustering

    pip install huggingface_hub llama-cpp-python
    pip install 'transformers[tf-cpu]'
