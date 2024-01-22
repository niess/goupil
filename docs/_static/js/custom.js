/* Add external link to commit hash. */
jQuery(document).ready(function(){
        var element = jQuery("span.commit").find("code").first();
        if (element.length) {
                const commit = element.text();
                element.replaceWith(
                    `<a href="https://github.com/niess/goupil/tree/${commit}">${commit}<a>`
                );
        }
});
