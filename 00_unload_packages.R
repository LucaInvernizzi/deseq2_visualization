lapply(names(sessionInfo()$otherPkgs), function(pkgs)
        detach(
                paste0('package:', pkgs),
                character.only = T,
                unload = T,
                force = T
        ))