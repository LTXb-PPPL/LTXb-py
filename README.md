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

- On github.com navigate to the LTXb-PPPL/LTXb-py repository
- Click the green `Code` dropdown button. The `HTTPS` option should be selected by default. Copy the repo link to the clipboard (should be https://github.com/LTXb-PPPL/LTXb-py.git).
- In PyCharm, using the menu bar, select VCS -> Git -> Clone
- Paste the repo link into the URL box. The Directory might default to some location on your local drive. It is recommended to use a mapped network drive. For example, I used Z:\PyCharmProjects\LTXb-py\ which is then located on NoMachine at /u/wcapecch/PyCharmProjects/LTXb-py/.
- Hit Clone. (Hitting the `Test` button first will check your connection to GitHub). This will create a local copy of the repo and allow you to make edits to the files. Any changes need to be committed and pushed before they'll be reflected in the main repo.
- *NOTE*: I believe you will also need to download and install git and I don't remember if that's something you need IT to do for you. It'd be great if edits were made to this by those who have to go through this for the first time.
- *Another NOTE*: PyCharm will ask you to login to GitHub at some point in the cloning process. I ran into issues using my username/password combo. If you experience a similar problem, try using the authentication token option instead.
