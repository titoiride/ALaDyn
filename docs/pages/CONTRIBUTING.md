title: Contributing guidelines

# Contributing guidelines

`ALaDyn` is an open source code, so you are very welcome to contribute!

## How to propose changes

According to the general rules of Open Source code development, it is not possible to push any changes directly in the main repository. However, any Pull Request is evaluated by the community and then merged into the **master** branch. If you are new to GitHub, in the following there are the required steps in order to contribute.

### 1. Fork the repository

The first required step is to fork [`ALaDyn`](https://github.com/ALaDyn/ALaDyn) repository. 

- Click on "Fork", on the top right of the page, to fork the main repository into your GitHub account. From this moment on, your forked repository share the same history with `ALaDyn`, but may have a different future development.
- Since the intent to contribute, we suggest (so, this is an _optional_ step) to set the main `ALaDyn` repository as a remote with
  ```
  git remote add upstream https://github.com/ALaDyn/ALaDyn
  ```
  This way, you can always keep your master branch _up-to-date_ with all the development, by doing
  ```
  git fetch upstream
  git merge upstream/master
  ```
### 2. Modify and commit the changes

In your repository, create a new branch (which we suggest you to call `dev/<yourname>/<nameofthebranch>`) and work on that. You should commit the changes and then push them on you online repository. The requested steps are therefore
```
git checkout dev/<yourname>/<nameofthebranch>

...modify the code...

git add <changedfile>
git commit -m 'Description of the changes'
git push
```

### 3. Open a Pull Request

After you pushed all the changes on your online repository, you can open a **pull request** by clicking on "pull request".
Your changes will be therefore available to the community that can eventually discuss on them and decide when to merge them into the main repository.

## Recommendations

Any pull request should

* Be as compact as possible, to facilitate the discussion and merging process. You can always open more than one pull request
* Contain correction and/or new fragments of code that can at least compile: you should try the compilation during your work. However, the Continuous Integration will check the pull request status as soon as you open it: any broken code cannot be merged.
* If possible, be already **squashed**. In particular, you should clean every confused or repeated commit, to reduce the number of them to the significant ones. For example (and this happens on daily basis, so don't worry!) an history such as
    
    ```
    <hash1> Fist try to fix bug
    <hash2> Second try to fix bug
    <hash3> Maybe this time works
    <hash4> Final fix
    ```

    should become

    ```
    <hash> Bug fixed
    ```