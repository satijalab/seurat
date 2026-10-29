// Navbar enhancements for pkgdown.
(() => {
  // Keep anchor and sidebar offsets in sync with the actual navbar height.
  const setupNavbarHeight = () => {
    const navbar = document.querySelector('.navbar');
    if (!navbar) return;

    const root = document.documentElement;
    const update = () => {
      const rem = parseFloat(getComputedStyle(root).fontSize);
      const height = navbar.getBoundingClientRect().height / rem;
      root.style.setProperty('--seurat-navbar-height', `${height}rem`);
    };
    new ResizeObserver(update).observe(navbar);
    update();
  };

  // Theme toggle: persist the selection and keep its icon and label in sync.
  const setupThemeToggle = () => {
    const link = document.getElementById('theme-toggle');
    if (!link) return;

    const root = document.documentElement;
    const isDark = () => root.getAttribute('data-bs-theme') === 'dark';
    const button = document.createElement('button');
    button.id = 'theme-toggle';
    button.type = 'button';
    button.className = 'nav-link';
    const icon = document.createElement('span');
    icon.setAttribute('aria-hidden', 'true');
    button.append(icon);
    link.replaceWith(button);

    const update = () => {
      const dark = isDark();
      const label = dark ? 'Switch to light mode' : 'Switch to dark mode';
      icon.className = dark ? 'fa fa-sun' : 'fa fa-moon';
      button.setAttribute('aria-label', label);
      button.title = label;
    };
    button.addEventListener('click', () => {
      const theme = isDark() ? 'light' : 'dark';
      localStorage.setItem('theme', theme);
      root.setAttribute('data-bs-theme', theme);
    });
    new MutationObserver(update).observe(root, {
      attributes: true, attributeFilter: ['data-bs-theme']
    });
    update();
  };

  // Vignette categories: native details controls with mouse-hover support.
  const setupVignetteMenu = () => {
    const menu = document.querySelector('[aria-labelledby="dropdown-vignettes"]');
    if (!menu || menu.querySelector('.vignette-submenu')) return;

    document.getElementById('dropdown-vignettes')?.setAttribute('data-bs-auto-close', 'outside');
    let group;
    for (const item of Array.from(menu.children)) {
      const heading = item.querySelector('.dropdown-header');
      if (heading) {
        const details = document.createElement('details');
        details.className = 'vignette-submenu';
        // Named details groups let the browser keep only one category open.
        details.name = 'vignettes';
        const summary = document.createElement('summary');
        summary.className = 'dropdown-item';
        summary.textContent = heading.textContent;
        group = document.createElement('ul');
        group.className = 'dropdown-menu';
        details.append(summary, group);
        item.replaceChildren(details);
        details.addEventListener('pointerenter', event => {
          if (event.pointerType !== 'mouse') return;
          details.open = true;
        });
        details.addEventListener('pointerleave', event => {
          if (event.pointerType === 'mouse') details.open = false;
        });
      } else if (item.querySelector('.dropdown-divider')) {
        item.remove();
      } else if (group) {
        group.append(item);
      }
    }
  };

  // Initialize after the navbar is available.
  const setup = () => {
    setupVignetteMenu();
    setupThemeToggle();
    setupNavbarHeight();
  };
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', setup, { once: true });
  } else {
    setup();
  }
})();
