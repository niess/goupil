jQuery(document).ready(function(){
        /* Add external link to commit hash. */
        var element = jQuery("li.wy-breadcrumbs-aside>a.fa.fa-github");
        if (element.length) {
                const base_url = element[0]
                        .href
                        .split("/blob/")[0];
                var element = jQuery("footer span.commit>code").first();
                if (element.length) {
                        const commit = element.text();
                        element.replaceWith(
                                `<a href="${base_url}/tree/${commit}">${commit}<a>`
                        );
                }
        }

        /* Add reference to Nord colour palette. */
        var element = jQuery("footer>a:contains('Read the Docs')")
                .first()
                .after(" and <a href='https://www.nordtheme.com/'>Nord</a> colours");
});
