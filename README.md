# LTXb-py
python routine repo

WELCOME!
This repo is intended to give a centralized, versioned location for python code that is widely used by the LTX-b group.
You should feel free to add code that you feel would be useful or edit code that is either broken or could be improved (discussing edits with the author is recommended).
Be sure to **comment** well and include mention of your authorship.

## Getting Started

### Access to the GitHub repo
You should have received an email inviting you to join the LTXb-PPPL GitHub Organization which led you here. If you want to invite new members contact Bill.

### Cloning the Repository
There's not much use to editing files online, so setting up a local clone is essential.

It is recommended that you clone to a mapped network drive instead of your local machine. While some code will work running from your local machine, many will not, for example routines making calls to MDSPlus or Discord bots.

#### PyCharm
As of setting this up, I use PyCharm 2019.2.2 (Community Edition) as my IDE and will walk through setting up a local clone of the LTXb-PPPL/LTXb-py repository. If people use other IDEs feel free to add specific instructions, but the steps here should give a general idea.

###### Software Requirements
- [PyCharm](https://www.jetbrains.com/pycharm/download/) (obviously)
- [Git](https://git-scm.com/downloads)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) (includes conda, Python, pip, and few others)

###### Cloning Procedure
- On github.com navigate to the LTXb-PPPL/LTXb-py repository
- Click the green `Code` dropdown button. The `HTTPS` option should be selected by default. Copy the repo link to the clipboard (should be https://github.com/LTXb-PPPL/LTXb-py.git).
- In PyCharm, using the menu bar, select VCS -> Git -> Clone
- Paste the repo link into the URL box. The Directory might default to some location on your local drive. It is recommended to use a mapped network drive. For example, I used Z:\PyCharmProjects\LTXb-py\ which is then located on NoMachine at /u/wcapecch/PyCharmProjects/LTXb-py/.
- Hit Clone. (Hitting the `Test` button first will check your connection to GitHub). This will create a local copy of the repo and allow you to make edits to the files. Any changes need to be committed and pushed before they'll be reflected in the main repo.
- *NOTE*: I believe you will also need to download and install git and I don't remember if that's something you need IT to do for you. It'd be great if edits were made to this by those who have to go through this for the first time.
- *Another NOTE*: PyCharm will ask you to login to GitHub at some point in the cloning process. I ran into issues using my username/password combo. If you experience a similar problem, try using the authentication token option instead.

###### Creating a Project Interpreter
Once you have the repository cloned PyCharm needs to establish an interpreter for the project. 

To create a new interpreter for this project (recommended)
- File -> Settings
- On the sidebar select `Project: LTXb-py` and click on `Project Interpreter`
- Click the gear icon next to the Interpreter box and select `Add`
- Select `Conda Environemnt` on the left, `New Environment` on the right, and select Python version 3.7
- Click `OK` to create a basic interpreter

The repo contains a requirements.txt file, so after assigning it an interpreter it will display a banner across the top of the code window if any dependencies are missing. Install any missing dependencies (this may take a while).

###### Updating requirements.txt
***NOTE***: This will document ***all*** the packages installed in your interpreter. To keep things clean, please create a unique interpreter for this project as detailed above and install ***only*** the dependencies necessary as listed in requirements.txt. If you then need new dependencies due to changes you've made to the repo you can update requirements.txt so that other users can easily update their interpreters.
To do this:
- open Anaconda prompt (anaconda into windows search bar)
- navigate to repo directory (Z:\PyCharmProjects\LTXb-py\ for example)
- `conda activate <interpreter name>`
- `pip freeze > requirements.txt`