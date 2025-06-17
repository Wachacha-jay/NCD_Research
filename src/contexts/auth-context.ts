
import { createContext } from 'react';
import { AuthContextType } from './auth-context-helpers';

export const AuthContext = createContext<AuthContextType | undefined>(undefined);
