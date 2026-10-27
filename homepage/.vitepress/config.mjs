import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Pseudovirus deep mutational scanning of Junin virus glycoprotein complex",
  description:
    "Experimental measurements of how mutations to the glycoprotein complex (GPC) from the Junin virus Romero strain affect cell entry and low pH stability.",
  base: "/JUNV-Romero-GPC-DMS/",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep/JUNV-Romero-GPC-DMS" }],
    footer: {
      message: 'See <a href="https://www.biorxiv.org/content/10.64898/2026.10.22.761119v1">Carr et al (2026)</a> for study details',
    },
  },
});
