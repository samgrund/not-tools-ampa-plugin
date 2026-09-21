# not-tools-ampa-plugin

This repository contains plugins for AMPA (Astro Multi-Purpose
Analyzer) specific to data from the Nordic Optical Telescope (NOT).

## Install

### Remote registry

In AMPA open **Plugins -> Plugin Browser -> Add Plugins -> Remote
Registry URLs** and add:

```
https://github.com/samgrund/not-tools-ampa-plugin/raw/main/plugins.json
```

Install the plugin from the browser list, then restart AMPA when
prompted.

### Local plugin directory

Clone the repository:

```console
$ git clone https://github.com/samgrund/not-tools-ampa-plugin.git
```

Add the cloned folder under **Plugins -> Plugin Browser -> Add
Plugins -> Local Plugin Directories**, then restart AMPA. Update with
`git pull` and a restart.
