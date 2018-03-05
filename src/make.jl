using Documenter, EAGOParametricInterval

makedocs(modules=[EAGOParametricInterval],
        doctest=true)

deploydocs(deps = Deps.pip("mkdocs", "python-markdown-math"),
    repo = "github.com/GITHUBNAME/GITHUBREPO.git",
    julia  = "0.6.0",
    osname = "windows")
