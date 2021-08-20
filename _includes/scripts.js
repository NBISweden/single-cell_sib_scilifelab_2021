var btns = document.getElementsByClassName("menuButton")
var url = window.location.href
for (var i = 0; i < btns.length; i++) {
  if (url  == btns[i].closest("a").href) {
    btns[i].classList.add('active')
  }
}
