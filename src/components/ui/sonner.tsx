
import { useTheme } from "next-themes"
import { Toaster as Sonner, toast } from "sonner"
import { toasterClassNames } from "./sonner-helpers"

type ToasterProps = React.ComponentProps<typeof Sonner>

const Toaster = ({ ...props }: ToasterProps) => {
  const { theme = "system" } = useTheme()

  return (
    <Sonner
      theme={theme as ToasterProps["theme"]}
      className="toaster group"
      toastOptions={{
        classNames: toasterClassNames,
      }}
      {...props}
    />
  )
}

export { Toaster, toast }
